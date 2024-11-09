#!/bin/bash
#SBATCH -J ragtag
#SBATCH -o ragtag.%j
#SBATCH -t 144:00:00
#SBATCH -N 1 -n 10
#SBATCH --mem=120G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

ragtag_scaffold.py /project/daane/hussain/final_project/raw_data/obeta_ncbi.fa /project/daane/hussain/final_project/hifiasm/obeta.asm.bp.p_ctg.fa -u
