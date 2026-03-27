# Intro

A Snakemake workflow for analyzing phage induction data. 

**NOTE**: geNomad, VirSorter2, Cenote-Taker3 and CheckV must be installed before running the workflow's identification step. The paths to these tools and their databases need to be updated to match your local setup by editing `phage-induction/config/config.yml`. Automatic installation will be added if the identification is kept.

## Set up environment

### Clone the repository 

```bash
git clone https://github.com/xwu35/phage-induction.git
```

### Install Snakemake

The workflow is built for Snakemake version 7. Version 8 and above introduce breaking changes and deprecations and have not been tested. It may not function correctly with newer versions. Please install Snakemake version 7 using the script below.

```bash
# option 1 using conda
conda env create -n snakemake -f phage-induction/snakemake_env.yml

# option 2 using mamba if it's installed
mamba env create -n snakemake -f phage-induction/snakemake_env.yml
```

If you have installed Snakemake version 7 on your system without using the commands above, make sure that `Biopython` is also installed. If it is not already installed, activate your snakemake environment and run `pip install biopython==1.86`.

### Download snakemake profile

The profile is required to run the workflow on HPC. Skip this step if you already have a SLURM profile in `~/.config/snakemake`.

```bash
# download the profile
git clone https://github.com/xwu35/slurm

# move the profile to the right directory
mv slurm ~/.config/snakemake 
```

### Export PATH

Add the phage-induction directory path to your environment variable so you can run `pinduction.py` without the full path.

```bash
echo 'export PATH="/your/path/to/phage-induction:$PATH"' >> ~/.bashrc
```

To make the changes take effect, you can either log out and log back in, or you can source the `.bashrc` file using the following command:

```bash
source ~/.bashrc
```

## Sample information table

The sample information table should look like this:
| sample                       | R1                                              | R2                                              | genome     | seq_name      | prophage_file                            |
|------------------------------|-------------------------------------------------|-------------------------------------------------|------------|---------------|------------------------------------------|
| 1041.01Kp_AMP_Digest_7_23_25 | 1041.01Kp_AMP_Digest_7_23_25_R1_catted.fastq.gz | 1041.01Kp_AMP_Digest_7_23_25_R2_catted.fastq.gz | 1041.01_KP | 1041.01_KP.fa | 1041.01_KP_final_coordinates_0-based.bed |
| 1041.01Kp_MMC_Digest_7_23_25 | 1041.01Kp_MMC_Digest_7_23_25_R1_catted.fastq.gz | 1041.01Kp_MMC_Digest_7_23_25_R2_catted.fastq.gz | 1041.01_KP | 1041.01_KP.fa | 1041.01_KP_final_coordinates_0-based.bed |
| 1144.01_E_AMP_Digest_8_28_25  | 1144.01_E_AMP_Digest_8_28_25_R1_catted.fastq.gz  | 1144.01_E_AMP_Digest_8_28_25_R2_catted.fastq.gz  | 1144.01_E | 1144.01_E.fa | 1144.01_E_final_coordinates_0-based.bed |
| 1144.01_E_AMP_Stock_8_12_25  | 1144.01_E_AMP_Stock_8_12_25_R1_catted.fastq.gz  | 1144.01_E_AMP_Stock_8_12_25_R2_catted.fastq.gz  | 1144.01_E | 1144.01_E.fa | 1144.01_E_final_coordinates_0-based.bed |

## Usage

The workflow supports two mapping software options (bowtie2 and minimap2), with bowtie2 used by default. Detailed usage information can be viewed using the -h or --help flags `pinduction.py -h`.

### Test run 

A dry-run can be performed to check which rules will be executed and which files will be produced. 

```bash
conda activate snakemake

pinduction.py --test  --dryrun
```

Do not run this on the login node. Submit it as an sbatch job. Check the HTCF usage guide here (https://github.com/xwu35/baldridge_lab/blob/main/HTCF.md). 

```bash
conda activate snakemake

pinduction.py --test  # An output directory named "test_output" will appear
```

### Example run

A dry-run can be performed by specifing the `--dryrun` flag

```bash
conda activate snakemake

pinduction.py \
    --sample_info path/to/sample/info/table \
    --reads_dir path/to/raw_reads_directory \
    --genome_dir path/to/genome_sequence_directory \
    --prophage_dir path/to/prophage_prediction_results \
    --output_dir output_dir \
    --step assemble 
```

### Specific steps

Specific steps can be run using the `--step` flag. 

- **fastqc**: QC on raw reads. It is not included in any other steps (default step)
- **trimming**: trim the reads and QC on trimmed reads 
- **mapping**: map trimmed reads back to the genome (includes trimming step)
- **assemble**: assemble the reads (includes trimming step)
- **alignment**: align the assembled contigs to the genome (includes trimming and assemble steps)
- **annotation**: functional annotation of the assembled contigs using pharokka (includes trimming and assemble steps)
- **all**: run all the steps above except for the fastqc step

The workflow runs fastqc by default.

## Output description

|             Directory          |                     Description                    |
|--------------------------------|----------------------------------------------------|
| `qc`                           | Includes results from fastqc and trimming          |
| `mapping`                      | Includes results from mapping reads to the genome  |
| `assembly`                     | Includes assembled contig sequences                | 
| `alignment`                    | Includes BLASTn results using the assembled contigs against the genome |
| `annotation`                   | Includes results from functional annotation on the assembled contigs |



