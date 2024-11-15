#!/bin/bash
#SBATCH -J chain_obeta_cgob
#SBATCH -o chain_obeta_cgob.%j
#SBATCH -t 120:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=10G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

module add Nextflow/21.10.6

python /project/daane/hussain/programs/make_lastz_chains/make_chains.py Cottoperca_gobio Opsanus_beta /project/daane/hussain/final_project/raw_data/cgob_dna_sm.fa /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi_sm.fa --project_dir chain_obeta_cgob --executor slurm --executor_queuesize 210 --seq1_chunk 50000000 --seq2_chunk 10000000

