# easy-NET

<p align="center">
  <img src="easy-NET.jpg" width="300"/>
</p>

**easy-NET** is a bash pipeline for constructing microbial co-occurrence networks using [SparCC](https://github.com/bio-developer/sparcc) (Sparse Correlations for Compositional data). It automates the full workflow — from ASV table filtering to correlation network generation and significance testing — making it accessible to researchers without extensive bioinformatics experience.

---

## Features

- **Flexible filtering**: filter ASV tables by absolute read count (`-n`), relative abundance (`-ab`), or let the script find the optimal cutoff automatically (`--auto`)
- **Automated network construction**: runs SparCC correlations for all `.txt` tables in a folder
- **Bootstrap analysis**: generates 100 permutations to assess correlation robustness
- **Significance testing**: computes two-sided pseudo p-values via permutation-based tests
- **Batch processing**: handles multiple tables in a single run, organizing outputs per table

---

## Requirements

- **Python 2.7** (required for SparCC compatibility)
- **SparCC** — install from [https://github.com/bio-developer/sparcc](https://github.com/bio-developer/sparcc)
- **conda** (recommended for environment management)

---

## Installation

Clone the repository:
```bash
git clone https://github.com/tpellegrinetti/easy-NET.git
cd easy-NET
```

Install and activate the conda environment:
```bash
bash install.sh
conda activate easy-NET
```

---

## Input Format

easy-NET expects tab-delimited `.txt` files with ASVs/OTUs in rows and samples in columns. Prepare one table per treatment or group of interest.

**Example:** if your experiment has three treatments, prepare:
```
data/
├── T1.txt
├── T2.txt
└── T3.txt
```

Each file should contain only the samples belonging to that treatment.

---

## Usage
```bash
bash easy_net.sh -d <folder> -s <sparcc_path> [--auto] [-n <min_reads>] [-ab <min_rel_abund_%>]
```

### Options

| Flag | Description |
|------|-------------|
| `-d <folder>` | Folder containing `.txt` ASV tables **(required)** |
| `-s <path>` | Path to the SparCC directory **(required)** |
| `--auto` | Automatically find the lowest cutoff keeping all tables < 1000 ASVs |
| `-n <int>` | Minimum absolute read count per ASV per sample |
| `-ab <float>` | Minimum relative abundance in % per sample (e.g. `0.1` = 0.1%) |
| `-h` | Show help message |

### Examples
```bash
# Automatic cutoff detection
bash easy_net.sh -d ./data -s /opt/SparCC --auto

# Filter by minimum read count
bash easy_net.sh -d ./data -s /opt/SparCC -n 100

# Filter by relative abundance
bash easy_net.sh -d ./data -s /opt/SparCC -ab 0.01

# Combine filters
bash easy_net.sh -d ./data -s /opt/SparCC -n 5 -ab 0.05
```

---

## Output

For each input table, easy-NET creates a `<table>_net/` directory containing:
```
T1_net/
├── cor_sparcc.out         # SparCC correlation matrix
├── pvals_two_sided.txt    # Two-sided pseudo p-values
├── perm/                  # 100 permuted tables
└── pvalues/               # SparCC correlations for each permutation
```

Filtered tables are saved in `<input_folder>/filtered/`.

---

## Workflow Overview
```
.txt tables → Filtering → SparCC correlations
                              ↓
                      Bootstrap permutations (n=100)
                              ↓
                      Pseudo p-value computation
                              ↓
                      Significant pair extraction
```

---

## Citation

If you use easy-NET in your research, please cite SparCC:

> Friedman J, Alm EJ (2012). Inferring Correlation Networks from Genomic Survey Data. *PLOS Computational Biology* 8(9): e1002687.

---

## License

MIT License. See `LICENSE` for details.
