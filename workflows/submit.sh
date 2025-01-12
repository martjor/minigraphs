#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=500M

module purge
module load Conda/3

conda init
conda activate minigraphs

snakemake --workflow-profile slurm --sdm conda --rerun-incomplete 