# easy-NET

easyNET is a user-friendly script designed to streamline the process of constructing microbial networks using SparCC (Sparse Correlations for Compositional data). This tool simplifies the workflow for analyzing compositional datasets and generating correlation-based microbial interaction networks, making it accessible for researchers without extensive programming expertise.

## Features:
### Data Filtering: Automatically filters input datasets based on absolute or relative abundance thresholds.
Automated Network Creation: Generates SparCC correlation networks with minimal user input.
### Bootstrap Analysis: Incorporates SparCC bootstrapping to assess the robustness of correlations.
Significance Testing: Computes p-values for network edges using permutation-based tests.
### Customizable Parameters: Allows users to set thresholds and specify directories for input and output data.

## Workflow:
### Input Processing: Reads tab-delimited .txt files containing microbial abundance data.
### Filtering: Removes rows based on user-defined abundance thresholds.

## Network Generation:
Calculates SparCC correlations.
Creates bootstrapped permutations for significance testing.
Generates p-value files to identify significant interactions.
### Output: Saves networks and significant edges in organized directories for easy interpretation.

## Why Use easyNET?
### Simplicity: Designed for researchers who need an efficient and straightforward way to create microbial networks.
### Automation: Handles complex SparCC workflows automatically, saving time and effort.
### Customization: Flexible options for input data filtering and network parameters.

# How to Use:
### Clone this repository:
git clone https://github.com/username/easyNET.git

### Install the script and activate conda:
conda env create -f environment.yml
conda activate easy-NET

### Run the script
./easyNET.sh -n <absolute_cutoff> -a <relative_cutoff> -d <data_directory>

# Dependencies:
Python 2.7 (for compatibility with SparCC)
SparCC: Ensure SparCC is installed and accessible in your system path. (https://github.com/bio-developer/sparcc)
With easyNET, you can focus on interpreting your results rather than troubleshooting complex pipelines. Simplify your microbial network analysis today!