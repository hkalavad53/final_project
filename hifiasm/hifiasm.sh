#!/bin/bash
#SBATCH -J hifiasm
#SBATCH -o hifiasm.%j
#SBATCH -t 144:00:00
#SBATCH -N 1 -n 10
#SBATCH --mem=120G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

/project/daane/hussain/programs/hifiasm/hifiasm -o obeta.asm -t 10 /project/daane/hussain/final_project/raw_data/obeta.fastq.gz

