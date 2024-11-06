# final_project

## Raw Data and Existing Data files

### Downloaded the chromosomal level assembly of Opsanus beta (Gulf Toadfish) from NCBI

```
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

