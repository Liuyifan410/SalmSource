# STMSource: Salmonella Typhimurium Host Origin Prediction

## Overview
STMSource is a bioinformatics tool for predicting the host origin (Chicken vs Pig) of *Salmonella* Typhimurium strains from genomic sequences. The tool implements a two-stage prediction pipeline. The first stage aims to distinguish isolates from chickens or swine from isolates of all other origins. The second stage then specifically discriminates between chicken and swine origins for isolates predicted to originate from one of these two hosts.

## Dependencies
All dependencies are managed through Conda:

**R Packages:**
- `tidymodels`
- `tidyverse` 
- `lightgbm`
- `bonsai`
- `seqinr`
- `ranger`

**Bioinformatics Tools:**
- `BLAST`
- `DIAMOND`
- `Snippy`
- `Prodigal`

## Installation

### Prerequisites
- Conda (Miniconda or Anaconda)

### Setup
```bash
git clone https://github.com/yourusername/STMSource.git
cd STMSource
conda env create -f environment.yml
conda activate STMSource
```

## Usage

### Basic Command Structure
```bash
bash SalmSource.sh -i <input_file> -p <output_prefix> -s <stage>
```

### Parameters

- `-i, --input`: Input genome file in FASTA format (required)
- `-p, --prefix`: Prefix for output files (required)  
- `-s, --stage`: Prediction stage - `one` or `two` (required)

## Output Files

### Stage One Outputs
- `*_groups.pep`: Protein sequences
- `*_groups.cds`: Coding sequences
- `*_groups_blast.out`: BLAST results
- `*_groups_host_prediction.txt`: Host prediction results
- `*_detected_features.txt`: Detected features

### Stage Two Outputs
- `*_target.pep`: Target protein sequences
- `*_target.cds`: Target coding sequences
- `*_blast.out`: BLAST results
- `*_host_prediction.txt`: Host prediction results
- `*_detected_features.txt`: Detected features
- `*_snippy/`: Snippy analysis directory

## License

This project is licensed under the terms of the LICENSE file included in the repository.

## Citation

If you use STMSource in your research, please cite this repository.



