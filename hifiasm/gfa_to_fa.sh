#!/bin/bash
#SBATCH -J gfa_to_fa
#SBATCH -o gfa_to_fa.%j
#SBATCH -t 144:00:00
#SBATCH -N 1 -n 10
#SBATCH --mem=120G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

awk '/^S/{print ">"$2;print $3}' /project/daane/hussain/final_project/hifiasm/obeta.asm.bp.p_ctg.gfa > obeta.asm.bp.p_ctg.fa
