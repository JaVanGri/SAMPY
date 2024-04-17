<img src="Logo_SAMPY.png" alt="Beschreibung" style="width: 100%;text-align:center">

# SAMPY: Python Implementation for Significance Analysis of Microarrays

SAMPY is a Python tool based on the Significance Analysis of Microarrays (SAM) method, designed to identify significant differences in measurements from multiple items across two different settings while taking the false discovery rate into account. This implementation follows the methodology outlined by Schwender in the reference "[Statistical Methods for Microarray Data Analysis: Methods and Protocols](https://link.springer.com/chapter/10.1007/3-540-28084-7_42)."


## How to use it?
To use SAMPI, follow these steps:

1. **Install the dependencies**:
   - Ensure you have Conda installed on your system.
   - All necessary dependencies are listed in `environment.yml`. Install them by running the following command in the terminal:
   ```
     conda env create -f environment.yaml
   ```

2. **Activate the Conda environment**:
   ```
     conda activate sampy_env
   ```
3. **Run the application**:
- Execute the application with the command:
  ```
  python outcyte_ups_v2.py /path/to/your/fasta/file.fasta
  ```

4. **Initial setup**:
- The first run might take a while since it needs to download a model from [Facebook Research's ESM](https://github.com/facebookresearch/esm).

5. **Access the results**:
- The result file will be stored at `/path/to/your/fasta/file.fasta.RESULT.csv`.


