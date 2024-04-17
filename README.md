<img src="Logo_SAMPY.png" alt="Beschreibung" style="width: 50%;text-align:center">

## What is SAMPY?
**OutCyte-UPS 2.0** is an innovative tool for predicting unconventional protein secretion based on the amino acid sequence of a protein.

## How to use it?
To use OutCyte-UPS 2.0, follow these steps:

1. **Install the dependencies**:
   - Ensure you have Conda installed on your system.
   - All necessary dependencies are listed in `environment.yml`. Install them by running the following command in the terminal:
   ```
     conda env create -f environment.yaml
   ```

2. **Activate the Conda environment**:
   ```
     conda activate oc2
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


