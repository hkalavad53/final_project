#!/bin/bash
#SBATCH -J repeatmask
#SBATCH -o repeatmask.o%j
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --mem=60G
#SBATCH -t 106:52:53

module load RepeatModeler/2.0.4-foss-2022a

module load RepeatMasker/4.1.5-foss-2022a

RepeatModeler -database obeta -threads 40 -LTRStruct

RepeatMasker -lib obeta-families.fa -pa 28 -e rmblast -xsmall -dir . -s obeta_ragtag_ncbi.fa

