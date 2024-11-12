#!/bin/bash
#SBATCH -J chain_<>_gacul
#SBATCH -o chain_<>_gacul.%j
#SBATCH -t 120:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=10G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

module add Nextflow/21.10.6

basedir=/project/daane/hussain/final_project/make_lastz_chains/obeta_ncbi

csh $basedir/chain_obeta_cgob/cleanUp.csh

gunzip $basedir/chain_obeta_cgob/Cottoperca_gobio.Opsanus_beta.allfilled.chain.gz

/project/daane/shared/TOGA/toga.py $basedir/chain_obeta_cgob/Cottoperca_gobio.Opsanus_beta.allfilled.chain /project/daane/hussain/final_project/TOGA/cgob.bed $basedir/chain_obeta_cgob/Cottoperca_gobio.2bit $basedir/chain_obeta_cgob/Opsanus_beta.2bit --pn toga_obeta_cgob -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt --nc /project/daane/shared/TOGA/nextflow_config_files --cesar_buckets 2,5,10,50

