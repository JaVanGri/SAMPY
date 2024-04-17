# %%

import argparse

import random
import math
import numpy as np
import scipy.interpolate as interpolate


def binomial_coefficient(n, k):
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))
def count_combinations(ones, zeros):
    return binomial_coefficient(ones + zeros, ones)

def calcSplitIdx(m:int,n:int,N:int=100,seed=1841):

    def generate_random_permutation(ones, zeros):
        permutation = [True] * ones + [False] * zeros
        random.shuffle(permutation)
        return tuple(permutation)  # Tupel sind hashable und können in einem Set gespeichert werden

    def get_k_unique_permutations(ones, zeros, K):
        permutations = set()
        while len(permutations) < K:
            perm = generate_random_permutation(ones, zeros)
            permutations.add(perm)
        return list(permutations)

    max_splits = count_combinations(m,n)
    splits = get_k_unique_permutations(m,n,min(N,max_splits))

    permutations=[]
    g1Index=np.array(list(range(m)))
    g2Index=np.array(list(range(n)))

    np.random.seed(seed)

    permCount=len(splits)
    permChoise= np.random.permutation(range(permCount))[0:min(N,permCount)]
    for i,split in enumerate(splits):
        if i in permChoise:
            comb1Split=np.array(split[0:m])
            comb2Split=np.array(split[m:m+n])

            comb_1 = g1Index[comb1Split]
            comb_2 = g2Index[comb2Split]
            permutations.append( [comb_1,comb_2])
    return permutations, len(permutations)

def d(state1: np.array, state2: np.array, s0: float):
    avg_state1 = state1.mean()
    avg_state2 = state2.mean()

    m = len(state1)
    n = len(state2)

    a = (1 / m + 1 / n) / (m + n - 2)

    s = np.sqrt(a * (np.sum((state1 - avg_state1) ** 2) + np.sum((state2 - avg_state2) ** 2)))

    d = (avg_state1 - avg_state2) / (s + s0)

    return d

def calculate_dvals(secretom, lysate, s0):
    d_list = []
    for i in range(len(secretom)):
        state1 = secretom[i]
        state2 = lysate[i]

        d_list.append(d(state1, state2, s0))
    return np.array(d_list)

def calcDvalDistWithPerm(group1, group2, s0=0.8, N=1000):
    """Calculate the distribution of d-values for permutations of the provided groups."""
    num_cols_group1, num_cols_group2 = group1.shape[1], group2.shape[1]
    perms, _ = calcSplitIdx(num_cols_group1, num_cols_group2, N=N)

    perm_d_values = []
    for perm in perms:
        # Generate permutations for group1 and group2
        perm_group1 = np.hstack((group1[:, perm[0]], group2[:, perm[1]]))
        remaining_indices_group1 = np.delete(np.arange(num_cols_group1), perm[0])
        remaining_indices_group2 = np.delete(np.arange(num_cols_group2), perm[1])
        perm_group2 = np.hstack((group1[:, remaining_indices_group1], group2[:, remaining_indices_group2]))

        # Calculate d-values for the current permutation
        d_values_perm = calculate_dvals(perm_group1, perm_group2, s0)
        perm_d_values.append(np.sort(d_values_perm))

    # Convert the list of permuted d-values into an array and compute the mean
    perm_d_values_array = np.array(perm_d_values).T
    expected_d_values = perm_d_values_array.mean(axis=1)

    return expected_d_values, perm_d_values_array


def calculate_pvals(dVals, dValsPermut):
    """Calculate p-values for observed d-values based on permutation d-values."""
    dValsPermut = np.abs(np.array(dValsPermut))
    m, B = dValsPermut.shape

    pVals = []
    for dVal in dVals:
        pVal = np.sum(np.abs(dValsPermut) >= np.abs(dVal)) / (m * B)
        pVals.append(pVal)
    return pVals

def calculate_PI0(pVals, lambdas=np.linspace(0, .95, 95)):
    """Estimate the proportion of true null hypotheses (pi0) from p-values."""
    m = len(pVals)
    piY = np.array([np.sum(pVals > l) / ((1 - l) * m) for l in lambdas])

    # Using UnivariateSpline for smoothing the pi0 estimates
    spline = interpolate.UnivariateSpline(lambdas, piY)

    # Evaluate the spline at lambda=0.5, limit the result between [1e-8, 1] to ensure valid probability values
    pi0_estimate = np.clip(np.mean(spline(0.5)), 1e-8, 1)

    return pi0_estimate



def calculate_FDR(observed_d_values, permuted_d_values, delta_threshold, pi0):
    """Calculate the False Discovery Rate (FDR) for a set of observed and permuted d-values."""
    permuted_d_values = np.array(permuted_d_values)
    m, B = permuted_d_values.shape
    sorted_d_values = np.sort(observed_d_values)
    expected_d_values = permuted_d_values.mean(axis=1)

    # Find the index of the smallest absolute deviation in expected values
    index_zero = np.argmin(np.abs(expected_d_values))

    # Determine upper cutoff
    upper_index = next((i for i in range(index_zero, m) if sorted_d_values[i] - expected_d_values[i] >= delta_threshold), m-1)
    upper_cutoff = sorted_d_values[upper_index] if upper_index < m else np.Inf

    # Determine lower cutoff
    lower_index = next((i for i in range(index_zero, -1, -1) if sorted_d_values[i] - expected_d_values[i] <= -delta_threshold), 0)
    lower_cutoff = sorted_d_values[lower_index] if lower_index >= 0 else -np.Inf

    # Counting the number of permuted values beyond the cutoffs
    count_above = np.sum(permuted_d_values >= upper_cutoff, axis=1)
    count_below = np.sum(permuted_d_values <= lower_cutoff, axis=1)
    numerator = np.sum(count_above + count_below)

    denominator = B * max(1, m + 1 - (upper_index - lower_index))

    return pi0 * numerator / denominator, lower_cutoff, upper_cutoff


def calcDeltaForFDR(dVals, dValsExp, dValsPermut, pi0, fdr_target=0.05, delta_min=0.1, fdr_tolerance=1e-2, max_iterations=500):
    """Calculate the delta value for False Discovery Rate adjustment using binary search."""
    delta_max = np.max(np.abs(np.sort(dVals) - np.sort(dValsExp)))
    left, right = delta_min, delta_max
    best_delta = (left + right) / 2
    count = 0

    current_fdr, _, _ = calculate_FDR(dVals, dValsPermut, best_delta, pi0)
    best_fdr = current_fdr

    while abs(current_fdr - fdr_target) > fdr_tolerance and count < max_iterations and right - left > fdr_tolerance:
        if current_fdr > fdr_target:
            left = best_delta
        else:
            right = best_delta

        best_delta = (left + right) / 2
        current_fdr, _, _ = calculate_FDR(dVals, dValsPermut, best_delta, pi0)
        if abs(current_fdr - fdr_target) < abs(best_fdr - fdr_target):
            best_fdr = current_fdr

        count += 1

    return best_delta


def merge_dataframes_on_id(df1, df2):
    """Merges two dataframes on the 'id' column."""
    return pd.merge(df1, df2, on='id')

def prepare_data(merged_data, df1, df2):
    """Prepares data arrays from merged data based on the original group dataframes."""
    group1_cols = [col for col in df1.columns if col != 'id']
    group2_cols = [col for col in df2.columns if col != 'id']
    return np.array(merged_data[group1_cols]), np.array(merged_data[group2_cols])

def append_analysis_results(data, d_values, d_values_exp, p_values):
    """Appends the results of the analysis to the dataframe."""
    data['dVals'] = d_values
    data['pVals'] = p_values
    data = data.sort_values(by='dVals')
    data['dEVals'] = d_values_exp
    return data[['id', 'dVals', 'dEVals', 'pVals']]


def decide_significance(data, delta):
    """Decides whether the difference between dVal and expexted dVal is significant."""
    data['significant_difference'] = [abs(d - de) > delta for d, de in zip(data['dVals'], data['dEVals'])]
    return data

def rename_columns(df, suffix):
    """Rename columns except 'id' with a specified suffix."""
    df = df.rename(columns={col: col + suffix for col in df.columns if col != 'id'})
    return df

def analyse(data_group1, data_group2, s0=0.8, fdr_target=0.05, fdr_tolerance=1e-2, max_delta_search_iterations=100, min_delta=0.01, max_permutations=1000):

    data_group1 = rename_columns(data_group1, '_1')
    data_group2 = rename_columns(data_group2, '_2')

    merged_data = merge_dataframes_on_id(data_group1, data_group2)
    group1_data, group2_data = prepare_data(merged_data, data_group1, data_group2)

    d_values = calculate_dvals(group1_data, group2_data, s0)
    d_values_exp, d_values_perm = calcDvalDistWithPerm(group1_data, group2_data, s0, N=max_permutations)

    p_values = calculate_pvals(d_values, d_values_perm)
    pi0 = calculate_PI0(p_values)
    delta = calcDeltaForFDR(d_values,
                            d_values_exp,
                            d_values_perm,
                            pi0,
                            fdr_target=fdr_target,
                            delta_min=min_delta,
                            max_iterations=max_delta_search_iterations,
                            fdr_tolerance=fdr_tolerance)
    fdr, _, _ = calculate_FDR(d_values, d_values_perm, delta, pi0)
    analysis_results = append_analysis_results(merged_data, d_values, d_values_exp, p_values)
    analysis_results = decide_significance(analysis_results,delta)

    return analysis_results, delta, fdr



import pandas as pd
import numpy as np
import sys

######################################################################
# Command Line Tool
######################################################################
def read_and_validate_csv(file_path):
    """ Read and validate the CSV data. """
    try:
        # Read the CSV file
        data = pd.read_csv(file_path)

        # Ensure there is an 'id' column
        if 'id' not in data.columns:
            raise ValueError("The CSV must contain an 'id' column.")

        # Check and convert all other columns to numeric types, replacing NaNs with 0
        for col in data.columns:
            if col != 'id':
                original_count = data[col].isna().sum()
                data[col] = pd.to_numeric(data[col], errors='coerce').fillna(0)
                if original_count > 0 or data[col].isna().sum() > 0:
                    print(f"Invalid values in '{col}' have been replaced with 0.")

        return data

    except Exception as e:
        print(f"Error reading the file {file_path}: {e}")
        sys.exit(1)

def main(args):
    data1 = read_and_validate_csv(args.file_path1)
    data2 = read_and_validate_csv(args.file_path2)

    analysis_results, delta, fdr = analyse(
        data1, data2,
        args.s0, args.fdr_target, args.fdr_tolerance,
        args.max_delta_search_iterations, args.min_delta, args.max_permutations
    )
    print("Analysis results:")

    analysis_results.to_csv(args.target_location + "SAMPY_result.csv",index=False)
    print(f"Found {np.sum(analysis_results[analysis_results['dVals'] - analysis_results['dEVals'] > 0]['significant_difference'])} positives and {np.sum(analysis_results[analysis_results['dVals'] - analysis_results['dEVals'] < 0]['significant_difference'])} negative differences with an FDR of {fdr:.2f}")
    print(f"Delta: {delta:.2f}")
    print(f"The result file is stored at {args.target_location + 'SAMPY_result.csv'}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two CSV files with statistical analysis.")
    parser.add_argument("file_path1", type=str, help="Path to the first CSV file")
    parser.add_argument("file_path2", type=str, help="Path to the second CSV file")
    parser.add_argument("--s0", type=float, default=0.8, help="Stability constant for d-value calculation")
    parser.add_argument("--fdr_target", type=float, default=0.05, help="Target value for the False Discovery Rate")
    parser.add_argument("--fdr_tolerance", type=float, default=0.01, help="Tolerance for the FDR calculation")
    parser.add_argument("--max_delta_search_iterations", type=int, default=100, help="Maximum number of iterations for delta search")
    parser.add_argument("--min_delta", type=float, default=0.01, help="Minimum delta value")
    parser.add_argument("--max_permutations", type=int, default=1000, help="Maximum number of permutations")
    parser.add_argument("--target_location", type=str, default="", help="Path where the result should be saved")

    args = parser.parse_args()
    main(args)
