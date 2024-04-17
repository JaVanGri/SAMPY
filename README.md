<img src="Logo_SAMPY.png" alt="Beschreibung" style="width: 100%;text-align:center">

# SAMPY: Python Implementation for Significance Analysis of Microarrays

SAMPY is a Python tool based on the Significance Analysis of Microarrays (SAM) method, designed to identify significant differences in measurements from multiple items across two different settings while taking the false discovery rate into account. This implementation follows the methodology outlined by Schwender in the reference "[Statistical Methods for Microarray Data Analysis: Methods and Protocols](https://link.springer.com/chapter/10.1007/3-540-28084-7_42)."


## How to use it?
To use SAMPI, follow these steps:

1. **Install the dependencies**:
   - Ensure you have Conda installed on your system.
   - All necessary dependencies are listed in `environment.yml`. Install them by running the following command in the terminal:
   ```
     conda env create -f environment.yml
   ```

2. **Activate the Conda environment**:
   ```
     conda activate SAMPY
   ```
3. **Run the application**:
   Use the command below to run SAMPY with the necessary arguments. Replace `path/to/file1.csv` and `path/to/file2.csv` with the paths to your input CSV files.

   ```bash
   python sampy.py path/to/file1.csv path/to/file2.csv --s0 0.8 --fdr_target 0.05 --fdr_tolerance 0.01 --max_delta_search_iterations 100 --min_delta 0.01 --max_permutations 1000
5.  arameter Explanation

Below is a detailed explanation of each parameter that can be used with the SAMPY tool:

- `file_path1`: **Required.** Path to the first input CSV file containing the dataset for the first group. This file should include an 'id' column and other columns with numeric values representing different measurements.

- `file_path2`: **Required.** Path to the second input CSV file containing the dataset for the second group. Similar to `file_path1`, this file must include an 'id' column and other numeric columns.

- `--s0` (default `0.8`): **Optional.** The stability constant used in the d-value calculation. A higher value increases the stability of the variance estimation, which can be useful in datasets with high variability.

- `--fdr_target` (default `0.05`): **Optional.** The target False Discovery Rate (FDR). This is the expected proportion of false positives among the declared significant results. Lowering this value tightens the criteria for significance, reducing the likelihood of false positives.

- `--fdr_tolerance` (default `0.01`): **Optional.** The tolerance for the convergence of the FDR calculation. This parameter determines how close the estimated FDR should be to the `fdr_target` during the iterative calculation process.

- `--max_delta_search_iterations` (default `100`): **Optional.** The maximum number of iterations allowed in the delta search process. This search is crucial for determining the threshold delta that achieves the target FDR. More iterations can provide a more precise result but at the cost of increased computational time.

- `--min_delta` (default `0.01`): **Optional.** The minimum possible value of delta used in the FDR calculations. This acts as a lower bound to the adjustments made for multiple testing.

- `--max_permutations` (default `1000`): **Optional.** The maximum number of permutations used to estimate the null distribution of the test statistic. Increasing this number can lead to more accurate and stable estimates of the p-values but requires more computation.

- `--target_location`: **Optional.** If specified, the path where the results of the analysis will be saved. If not provided, results will be output to the console.

Each parameter is designed to offer flexibility in adjusting the analysis to fit various datasets and experimental designs, providing robustness and adaptability in the application of the SAM method.




