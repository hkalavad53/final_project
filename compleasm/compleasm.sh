#!/bin/bash
#SBATCH -J compleasm
#SBATCH -o compleasm.%j
#SBATCH -t 144:00:00
#SBATCH -N 1 -n 10
#SBATCH --mem=25G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

compleasm run -a /project/daane/hussain/final_project/raw_data/obeta_ncbi.fa -o compleasm_obeta_ncbi -l actinopterygii_odb10

compleasm run -a /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi.fa -o compleasm_obeta_ragtag_ncbi -l actinopterygii_odb10

compleasm run -a /project/daane/hussain/final_project/raw_data/obeta_ragtag_tama.fa -o compleasm_obeta_ragtag_tama -l actinopterygii_odb10

