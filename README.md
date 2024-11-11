# final_project (on Sabine cluster)

## Raw Data and Existing Data files (not committed to repository due to huge file sizes)

### Downloaded the chromosomal level assembly of Opsanus beta (Gulf Toadfish) from NCBI

```
mkdir /project/daane/hussain/final_project

cd /project/daane/hussain/final_project

mkdir raw_data

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/802/735/GCA_040802735.1_UM_OBeta_2.1/GCA_040802735.1_UM_OBeta_2.1_genomic.fna.gz

gunzip https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/802/735/GCA_040802735.1_UM_OBeta_2.1/GCA_040802735.1_UM_OBeta_2.1_genomic.fna.gz

mv GCA_040802735.1_UM_OBeta_2.1_genomic.fna obeta_ncbi.fna

``` 

### Copied existing data files for O. beta hifi reads

``` 
cp /project/daane/hussain/assemblies/o_beta/hifi_reads/m84100_231231_051457_s3.hifi_reads.bc1002.bam /project/daane/hussain/final_project/raw_data/.

cp /project/daane/hussain/assemblies/o_beta/hifi_reads/m84100_231231_051457_s3.hifi_reads.bc1002.bam.pbi /project/daane/hussain/final_project/raw_data/.
```

These reads were scaffolded using ragtag's reference-guided scaffolding method using Thallasophryne amazonica (a close relative) as the reference. Now for this project I will use ragtag to scaffold these reads using the O. beta genome available on NCBI to check for differences in the genome quality.

```
cp /project/daane/hussain/assemblies/o_beta/scaffold/ragtag_output/ragtag.scaffold.fasta /project/daane/hussain/final_project/raw_data/obeta_ragtag_tamaz.fa
```

## Scaffold Hifi Reads to Opsanus beta genome obtained from NCBI

### Convert .bam to .fastq

```
sbatch bam2fastq.sh
```
``` 
#!/bin/bash
#SBATCH -J bam2fastq
#SBATCH -o bam2fastq.%j
#SBATCH -t 14:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=10G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

bam2fastq /project/daane/hussain/final_project/raw_data/m84100_231231_051457_s3.hifi_reads.bc1002.bam -o obeta
```

### Run Hifiasm to assemble the hifi reads

```
cd /project/daane/hussain/final_project
mkdir hifiasm
cd hifiasm
sbatch hifiasm.sh
```
```
#!/bin/bash
#SBATCH -J hifiasm
#SBATCH -o hifiasm.%j
#SBATCH -t 144:00:00
#SBATCH -N 1 -n 10
#SBATCH --mem=120G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

/project/daane/hussain/programs/hifiasm/hifiasm -o obeta.asm -t 10 /project/daane/hussain/final_project/raw_data/obeta.fastq.gz
```

### Run Ragtag to scaffold the assembly on Chromosomal-level O. beta genome obtained from NCBI

```
cd /project/daane/hussain/final_project
mkdir ragtag
cd ragtag
sbatch ragtag.sh
```
```
#!/bin/bash
#SBATCH -J ragtag
#SBATCH -o ragtag.%j
#SBATCH -t 144:00:00
#SBATCH -N 1 -n 10
#SBATCH --mem=120G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

ragtag_scaffold.py /project/daane/hussain/final_project/raw_data/obeta_ncbi.fa /project/daane/hussain/final_project/hifiasm/obeta.asm.bp.p_ctg.fa -u
```

### Run Compleasm on all assemblies to see differences in quality

```
cd /project/daane/hussain/final_project
mkdir compleasm
cd compleasm
wget https://busco-archive.ezlab.org/data/lineages/actinopterygii_odb10.2024-01-08.tar.gz
tar -xzf actinopterygii_odb10.2024-01-08.tar.gz
rm -r actinopterygii_odb10.2024-01-08.tar.gz
sbatch compleasm.sh
```
```
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
```

## Create Pairwise Sequence Alignment with Cottoperca gobio genome using make_lastz_chains pipeline from Michael Hiller Lab

### Check for softmask percentage

```
cd /project/daane/hussain/final_project
mkdir repeatmask
cd repeatmask
sbatch softmask_percentage.sh
```
```
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
```

The genome obtained from NCBI is already softmasked and the ragtag assembly using T. amazonica as reference has previously been softmasked and ran through the make_lastz_chains pipeline. Here, I will use windowmasker to softmask repeat sequences in ragtag assembly which used NCBI's O. beta genome as reference.

```
sbatch windowmasker.sh
```
```
#
```
