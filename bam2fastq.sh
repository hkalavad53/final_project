#!/bin/bash
#SBATCH -J bam2fastq
#SBATCH -o bam2fastq.%j
#SBATCH -t 14:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=10G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

bam2fastq raw_data/m84100_231231_051457_s3.hifi_reads.bc1002.bam -o obeta

