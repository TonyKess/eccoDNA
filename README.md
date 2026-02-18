# eccoDNA: eDNA Metabarcoding Workflow
![eccoDNA](https://github.com/user-attachments/assets/90abf742-23e1-45a3-a669-2621e6d11bdb)

A Nextflow DSL2 pipeline for processing environmental DNA samples through quality control, denoising, chimera removal, and taxonomic assignment.

## Installation

```bash
git clone https://github.com/TonyKess/eccoDNA.git
cd eccoDNA
./download_apptainer_images.sh
```

This downloads all required Apptainer containers to `./containers`.

## Setup
This workflow is meant to be run with independent samples and markers, matching the data structure of sequencing output for multi-marker eDNA projects. The directory structure is expected to look like the example below, where each directory
is a separate sequencing run of a unique eDNA primer.

```
project/
├── FISHE/
│   ├── samplesheet.csv
│   └── nextflow.FISHE.config
├── MIFISHU/
│   ├── samplesheet.csv
│   └── nextflow.MIFISHU.config
├── MCINNES16S/
│   ├── samplesheet.csv
│   └── nextflow.MCINNES16S.config

```
### 1. Create Sample Sheet
First, copy the helper scripts to your working directory:
```bash
cp creates_samplesheet.sh project/
```

The `create_samplesheet.sh` script generates a CSV from paired-end FASTQ files:

```bash
./create_samplesheet.sh /path/to/fastq/files
```

Or if in the project and marker directories already:

```bash
./create_samplesheet.sh .
```

Output: `samplesheet.csv` with columns: `sample,read1,read2`

**Process multiple directories:**

```bash
#!/bin/bash
for DIR in FISHE MIFISHU MCINNES16S ; do
    echo "Processing $DIR..."
    ./create_samplesheet.sh "$DIR"
    mv samplesheet.csv  $DIR ;
done
```

### 2. Configure Per Project

In your project directory, copy a nextflow config to each directory (you can rename them if you'd like):

```bash
mkdir project/configs
cp nextflow.generic.config project/FISHE/nextflow.config
cp nextflow.generic.config project/MIFISHU/nextflow.config
cp nextflow.generic.config project/MCINNES16S/nextflow.config
```

Edit each `project/configs/nextflow.MARKER.config` and change:
- `project_name` - descriptive name
- `marker` - marker name (FISHE, MIFISHU, etc.)
- `container_base` - path to containers directory
- `blast.database` - path to BLAST database
- `blast.db` - BLAST database name
- `clusterOptions` - SLURM account, partition, comment

### 3. Configure Job Submission
Copy the helper scripts to your project directory
```bash
cp submit.eccodna.customconfig.sh project/FISHE
cp submit.eccodna.customconfig.sh project/MIFISHU
cp submit.eccodna.customconfig.sh project/MCINNES16S
```

Edit `submit.eccodna.customconfig.sh`:
- `--account`, `--partition`, `--comment` - SLURM options
- `NXF_TEMP`, `APPTAINER_TMPDIR`, `TMPDIR` - temp directory paths
- `PIPELINE_SCRIPT` - path to main pipeline directory

### 4. Submit Jobs

For any marker, in the project/marker directory, for HPC SLURM submission do:

```bash
sbatch submit.eccodna.customconfig.sh
```
If custom configs have been renamed from nextflow.config, they can be specified directly in job submission with
```bash
sbatch --export=ALL,CONFIG_FILE=nextflow.FISHE.config submit.eccodna.customconfig.sh
'''

## Outputs

Results are organized in `${params.outdir}` with subdirectories:
- `merged/` - read merging logs
- `filtered/` - quality-filtered sequences
- `denoised/` - ASVs after denoising
- `final/` - ASV table, BIOM file, BLAST results, taxonomy assignments

## Options

Add to command line or edit config:

```bash
--skip_blast          # Skip taxonomic assignment
--skip_sss_filter     # Skip SSS filtering
--sss.min_threshold   # SSS threshold (default: 99.0)
```

## Development

Built as a Nextflow implementation of the eDNA analysis framework described in 
[Crowley et al. 2024.](https://doi.org/10.1002/edn3.517). Code structure was templated 
using Claude 4.5 Sonnet/Haiku, Gemini 3, ChatGPT 5.0.
