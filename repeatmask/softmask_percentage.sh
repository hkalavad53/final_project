#!/bin/bash
#SBATCH -J softmask_percentage
#SBATCH -o softmask_percentage.%j
#SBATCH -t 144:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=25G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

python ./softmask_percentage.py /project/daane/hussain/final_project/raw_data/obeta_ncbi.fa

python ./softmask_percentage.py /project/daane/hussain/final_project/raw_data/obeta_ragtag_tama.fa

python ./softmask_percentage.py /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi.fa

