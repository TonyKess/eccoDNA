# eccoDNA: eDNA Metabarcoding Workflow
![eccoDNA](![eccoDNA](https://github.com/user-attachments/assets/b42f16e7-d59f-40ed-b919-88388b3eedb2))

A Nextflow DSL2 pipeline for processing environmental DNA samples through quality control, denoising, chimera removal, and taxonomic assignment.

## Installation

```bash
git clone <eccoDNA-repo>
cd eccoDNA
./download_apptainer_images.sh
```

This downloads all required Apptainer containers to `./containers`.

## Setup

### 1. Create Sample Sheet

The `create_samplesheet.sh` script generates a CSV from paired-end FASTQ files:

```bash
./create_samplesheet.sh /path/to/fastq/files
```

Output: `samplesheet.csv` with columns: `sample,read1,read2`

**Process multiple directories:**

```bash
#!/bin/bash
for DIR in /data/project1 /data/project2 /data/project3; do
    echo "Processing $DIR..."
    ./create_samplesheet.sh "$DIR"
    mv samplesheet.csv samplesheet_$(basename $DIR).csv
done
```

### 2. Configure Per Project

Create a project directory with marker-specific subdirectories:

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
└── configs/
    ├── nextflow.FISHE.config
    ├── nextflow.MIFISHU.config
    └── nextflow.MCINNES16S.config
```

Create configs and populate them:

```bash
mkdir -p project/{FISHE,MIFISHU,MCINNES16S}
mkdir -p project/configs
cp nextflow.generic.config project/configs/nextflow.FISHE.config
cp nextflow.generic.config project/configs/nextflow.MIFISHU.config
cp nextflow.generic.config project/configs/nextflow.MCINNES16S.config
```

Edit each `project/configs/nextflow.MARKER.config` and change:
- `workDir` - set to project-specific path (e.g., `/path/to/project/FISHE`)
- `project_name` - descriptive name
- `marker` - marker name (FISHE, MIFISHU, etc.)
- `container_base` - path to containers directory
- `blast.database` - path to BLAST database
- `blast.db` - BLAST database name
- `clusterOptions` - SLURM account, partition, comment

### 3. Configure Job Submission

Edit `submit.eccodna.customconfig.sh`:
- `--account`, `--partition`, `--comment` - SLURM options
- `NXF_TEMP`, `APPTAINER_TMPDIR`, `TMPDIR` - temp directory paths
- `PIPELINE_SCRIPT` - path to main pipeline directory

### 4. Submit Jobs

Run all configured markers at once:

```bash
./joblauncher.eccodna.sh
```

The launcher searches `project/configs/` for all `.config` files and submits a job for each using `submit.eccodna.customconfig.sh`.

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
