![Compleasm](https://github.com/user-attachments/assets/0218380e-0451-43f5-8e20-867feaefd51c)# final_project (on Sabine cluster)

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

Compleasm results suggest that the genome completeness of two ragtag assemblies were not very different from each other. The NCBI reference genome has more single-copy complete genes compared to ragtag assemblies.
![Uploading Compleasm<svg width="1280" height="720" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" overflow="hidden"><g><rect x="0" y="0" width="1280" height="720" fill="#FFFFFF"/><g><path d="M303.429 383.824C374.982 383.824 432.988 441.83 432.988 513.384 432.988 584.938 374.982 642.943 303.429 642.943 254.044 642.943 208.954 614.869 187.164 570.552L303.429 513.384Z" fill="#0E2841" fill-rule="evenodd"/></g><g><path d="M187.164 570.552C155.591 506.34 182.049 428.692 246.26 397.119 261.02 389.861 277.026 385.48 293.425 384.211L303.429 513.384Z" fill="#E97132" fill-rule="evenodd"/></g><g><path d="M293.425 384.211C294.155 384.154 294.886 384.104 295.618 384.06L303.429 513.384Z" fill="#47D45A" fill-rule="evenodd"/></g><g><path d="M295.618 384.06 303.429 513.384Z" fill="#0F9ED5" fill-rule="evenodd"/></g><g><path d="M295.618 384.06C298.218 383.903 300.823 383.824 303.428 383.824L303.429 513.384Z" fill="#A02B93" fill-rule="evenodd"/></g><g><path d="M294.5 384.5 289.5 370.5 283.5 370.5" stroke="#A6A6A6" stroke-linejoin="round" stroke-miterlimit="10" fill="none"/></g><g><path d="M299.5 383.5 302.5 369.5" stroke="#A6A6A6" stroke-linejoin="round" stroke-miterlimit="10" fill="none"/></g><g><text fill="#404040" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(423.596 601)">67.72</text></g><g><text fill="#404040" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(128.902 445)">31.04</text></g><g><text fill="#404040" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(241.64 371)">0.27</text></g><g><text fill="#404040" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(304.161 369)">0.96</text></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="27" transform="translate(34.9219 339)">O. beta Ragtag Assembly using NCBI Assembly</text></g><g><rect x="210" y="689" width="11" height="12" fill="#0E2841"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(226.781 701)">S</text></g><g><rect x="251" y="689" width="12" height="12" fill="#E97132"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(268.214 701)">D</text></g><g><rect x="295" y="689" width="12" height="12" fill="#47D45A"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(312.206 701)">F</text></g><g><rect x="336" y="689" width="11" height="12" fill="#0F9ED5"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(352.751 701)">I</text></g><g><rect x="371" y="689" width="11" height="12" fill="#A02B93"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(387.664 701)">M</text></g><g><path d="M940 383.824C1011.55 383.824 1069.56 441.83 1069.56 513.384 1069.56 584.938 1011.55 642.943 940 642.943 891.368 642.943 846.832 615.709 824.672 572.42L940 513.384Z" fill="#0E2841" fill-rule="evenodd"/></g><g><path d="M824.672 572.42C792.068 508.726 817.27 430.661 880.964 398.056 895.743 390.491 911.839 385.836 928.375 384.346L940 513.384Z" fill="#E97132" fill-rule="evenodd"/></g><g><path d="M928.375 384.346C929.483 384.247 930.592 384.161 931.702 384.09L940 513.384Z" fill="#47D45A" fill-rule="evenodd"/></g><g><path d="M931.702 384.09 940 513.384Z" fill="#0F9ED5" fill-rule="evenodd"/></g><g><path d="M931.702 384.09C934.465 383.913 937.232 383.824 940 383.824L940 513.384Z" fill="#A02B93" fill-rule="evenodd"/></g><g><path d="M930.5 384.5 918.5 371.5 912.5 371.5" stroke="#A6A6A6" stroke-linejoin="round" stroke-miterlimit="10" fill="none"/></g><g><path d="M935.5 383.5 943.5 368.5 949.5 368.5" stroke="#A6A6A6" stroke-linejoin="round" stroke-miterlimit="10" fill="none"/></g><g><text fill="#404040" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(1060.74 600)">67.47</text></g><g><text fill="#404040" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(776.025 447)">31.1</text></g><g><text fill="#404040" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(870.179 371)">0.41</text></g><g><text fill="#404040" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(951.444 369)">1.02</text></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="27" transform="translate(626.79 339)">O. beta Ragtag Assembly using</text></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="27" transform="translate(985.303 339)">T. amazonica</text></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="27" transform="translate(1140.47 339)">Assembly</text></g><g><rect x="846" y="689" width="12" height="12" fill="#0E2841"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(863.353 701)">S</text></g><g><rect x="888" y="689" width="11" height="12" fill="#E97132"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(904.785 701)">D</text></g><g><rect x="932" y="689" width="11" height="12" fill="#47D45A"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(948.777 701)">F</text></g><g><rect x="972" y="689" width="12" height="12" fill="#0F9ED5"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(989.323 701)">I</text></g><g><rect x="1007" y="689" width="12" height="12" fill="#A02B93"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(1024.24 701)">M</text></g><g><path d="M0 0 170 75" stroke="#A6A6A6" stroke-linejoin="round" stroke-miterlimit="10" fill="none" fill-rule="evenodd" transform="matrix(1 0 0 -1 597.5 173.5)"/></g><g><path d="M597.5 184.5 767.5 259.5" stroke="#A6A6A6" stroke-linejoin="round" stroke-miterlimit="10" fill="none" fill-rule="evenodd"/></g><g><path d="M597.777 184.871C594.541 249.612 539.434 299.471 474.693 296.234 409.952 292.998 360.093 237.891 363.329 173.15 366.566 108.409 421.673 58.5497 486.414 61.7863 546.627 64.7966 594.767 112.936 597.777 173.15L480.553 179.01Z" fill="#0E2841" fill-rule="evenodd"/></g><g><path d="M891.24 183.406C888.812 231.961 847.482 269.356 798.927 266.928 764.077 265.186 733.547 243.024 721.092 210.43L803.322 179.01Z" fill="#E97132" fill-rule="evenodd"/></g><g><path d="M721.092 210.43C709.645 180.471 715.397 146.657 736.107 122.168L803.322 179.01Z" fill="#47D45A" fill-rule="evenodd"/></g><g><path d="M736.107 122.168C738.354 119.511 740.756 116.99 743.302 114.617L803.322 179.01Z" fill="#0F9ED5" fill-rule="evenodd"/></g><g><path d="M743.302 114.618C778.865 81.4692 834.567 83.427 867.715 118.99 883.942 136.399 892.428 159.637 891.24 183.406L803.322 179.01Z" fill="#A02B93" fill-rule="evenodd"/></g><g><path d="M597.777 173.15C597.973 177.055 597.973 180.966 597.777 184.871L480.553 179.01Z" fill="#196B24" fill-rule="evenodd"/></g><g><path d="M573.5 174.5 574.5 156.5" stroke="#A6A6A6" stroke-linejoin="round" stroke-miterlimit="10" fill="none"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(300.958 186)">98.41</text></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(809.426 291)">0.69</text></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(667.377 168)">0.27</text></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(691.599 109)">0.03</text></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(846.691 88)">0.6</text></g><g><text fill="#FFFFFF" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(532.156 156)">1.59</text></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="27" transform="translate(376.45 39)">O. beta genome obtained from NCBI</text></g><g><rect x="919" y="108" width="12" height="11" fill="#0E2841"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(936.193 119)">S</text></g><g><rect x="919" y="141" width="12" height="12" fill="#E97132"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(936.193 153)">D</text></g><g><rect x="919" y="175" width="12" height="12" fill="#47D45A"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(936.193 187)">F</text></g><g><rect x="919" y="209" width="12" height="12" fill="#0F9ED5"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(936.193 221)">I</text></g><g><rect x="919" y="243" width="12" height="12" fill="#A02B93"/></g><g><text fill="#595959" font-family="Aptos,Aptos_MSFontService,sans-serif" font-weight="400" font-size="21" transform="translate(936.193 255)">M</text></g></g></svg>.svg…]()



## Create Pairwise Sequence Alignment with Cottoperca gobio genome using make_lastz_chains pipeline from Michael Hiller Lab

### Check for softmask percentage

```
cd /project/daane/hussain/final_project
mkdir repeatmask
cd repeatmask
cp ~/Python/softmask_percentages.py .
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

#### Opsanus beta genome scaffolded using T. amazonica as reference (using Ragtag)

The ragtag assembly using T. amazonica as reference has previously been softmasked and ran through the make_lastz_chains pipeline. 


#### Opsanus beta genome obtained from NCBI:

The genome obtained from NCBI is already softmasked.


#### Opsanus beta genome scaffolded using NCBI's O.beta genome as reference (using Ragtag)

Here, I will use Repeatmasker to softmask repeat sequences in ragtag assembly which used NCBI's O. beta genome as reference.

```
mkdir obeta_ragatg_ncbi
cd obeta_ragatg_ncbi
cp /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi.fa .
module load RepeatModeler/2.0.4-foss-2022a
BuildDatabase --name obeta obeta_ragtag_ncbi.fa
sbatch repeatmask.sh
```
```
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
```

```
cp /project/daane/hussain/final_project/repeatmask/obeta_ragatg_ncbi/obeta_ragtag_ncbi.fa.masked /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi_sm.fa
```

### Run make_lastz_chains pipeline by Michael Hiler Lab to generate pairwise alignments between O. beta and Cottoperca gobio

```
cd /project/daane/hussain/final_project/
mkdir make_lastz_chains
cd make_lastz_chains/
mkdir obeta_ncbi
mkdir obeta_ragtag_tama
mkdir obeta_ragtag_ncbi
```

The fasta headers were cleaned prior to running make_lastz_chains as periods and spaces are a big NO-NO in the fasta headers for these programs to process them.

```
python /project/daane/hussain/final_roject/make_lastz_chains/fasta_rename_remove_period.py /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi_sm.fa /project/daane/hussain/final_project/raw_data/out.fa
```
```
import sys

input_file = sys.argv[1]
output_file= sys.argv[2]

with open(input_file, "r") as input_f, open(output_file, "w") as output_f:
    for line in input_f:
        if line.startswith(">"):
            header_components = line.strip().split(".")
            header_components = header_components[:1]
            modified_header =" ".join(header_components)
            output_f.write(modified_header + "\n")
        else:
            output_f.write(line)
```

```
mv /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi_sm.fa /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi_sm.fa.orig
mv /project/daane/hussain/final_project/raw_data/out.fa /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi_sm.fa
```


#### Opsanus beta genome scaffolded using T. amazonica as reference (using Ragtag)

make_lastz_chains has previously been run on this genome so I am copying the files here but not staging them on GitHub due to huge file size.

```
cd /project/daane/hussain/final_project/obeta_ragtag_tama/
cp -r /project/daane/hussain/repeatmask/batrachoidiformes/opsanus_beta/chain_obet_cgob/ .
```

#### Opsanus beta genome obtained from NCBI:

```
cd /project/daane/hussain/final_project/make_lastz_chains/obeta_ncbi/
touch chain_obet_cgob.sh
cp /project/daane/hussain/repeatmask/notothenioids/cottoperca_gobio/cgob_dna_sm.fa /project/daane/hussain/final_project/raw_data/.
# make_lastz_chains pipeline doesn't like periods or just long headers in general in fasta files
python fasta_rename_remove_period.py /project/daane/hussain/final_project/raw_data/obeta_ncbi.fa obeta_ncbi_sm.fa
sbatch chain_obet_cgob.sh
```
```
#!/bin/bash
#SBATCH -J chain_obeta_cgob
#SBATCH -o chain_obeta_cgob.%j
#SBATCH -t 120:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=10G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

module add Nextflow/21.10.6

python /project/daane/hussain/programs/make_lastz_chains/make_chains.py Cottoperca_gobio Opsanus_beta /project/daane/hussain/final_project/raw_data/cgob_dna_sm.fa /project/daane/hussain/final_project/make_lastz_chains/obeta_ncbi/obeta_ncbi_sm.fa --project_dir chain_obeta_cgob --executor slurm --executor_queuesize 210 --seq1_chunk 50000000 --seq2_chunk 10000000
```
### Opsanus beta genome scaffolded using NCBI's O.beta genome as reference (using Ragtag)

```
cd /project/daane/hussain/final_project/make_lastz_chains/obeta_ragtag_ncbi/
touch chain_obet_cgob.sh
python fasta_rename_remove_period.py /project/daane/hussain/final_project/raw_data/obeta_ragtag_ncbi_sm.fa obeta_ragtag_ncbi_sm.fa
sbatch chain_obet_cgob.sh
```
```
#!/bin/bash
#SBATCH -J chain_obeta_cgob
#SBATCH -o chain_obeta_cgob.%j
#SBATCH -t 120:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=10G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

module add Nextflow/21.10.6

python /project/daane/hussain/programs/make_lastz_chains/make_chains.py Cottoperca_gobio Opsanus_beta /project/daane/hussain/final_project/raw_data/cgob_dna_sm.fa /project/daane/hussain/final_project/make_lastz_chains/obeta_ncbi/obeta_ncbi_sm.fa --project_dir chain_obeta_cgob --executor slurm --executor_queuesize 210 --seq1_chunk 50000000 --seq2_chunk 10000000
```

## Identify Orthologous Genes between O. beta and C. gobio using TOGA (by Michael Hiller Lab)

### Opsanus beta genome scaffolded using T. amazonica as reference (using Ragtag)

TOGA has previously been run on this genome so I am copying the files here and staging some result files on github.

```
cd /project/daane/hussain/final_project
mkdir TOGA
cd TOGA
mkdir obeta_ragtag_tama
cd obeta_ragtag_tama
cp -r /project/daane/hussain/repeatmask/batrachoidiformes/opsanus_beta/toga_obeta_cgob .
```


### Opsanus beta genome scaffolded using NCBI's O.beta genome as reference (using Ragtag)

```
cd /project/daane/hussain/final_project/TOGA
cp /project/daane/hussain/repeatmask/notothenioids/cottoperca_gobio/cgob.bed .
cp /project/daane/hussain/repeatmask/notothenioids/cottoperca_gobio/cgob_isoforms.txt .
mkdir obeta_ncbi
cd obeta_ncbi
touch toga_obeta_cgob.sh
sbatch toga_obeta_cgob.sh
```
```
#!/bin/bash
#SBATCH -J chain_obeta_cgob
#SBATCH -o chain_obeta_cgob.%j
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
```


### Opsanus beta genome obtained from NCBI

```
cd /project/daane/hussain/final_project/TOGA
mkdir obeta_ragtag_ncbi
cd obeta__ragtag_ncbi
touch toga_obeta_cgob.sh
sbatch toga_obeta_cgob.sh
```
```
#!/bin/bash
#SBATCH -J chain_obeta_cgob
#SBATCH -o chain_obeta_cgob.%j
#SBATCH -t 120:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=10G
#SBATCH --mail-user=hskalavad@gmail.com
#SBATCH --mail-type=ALL

module add Nextflow/21.10.6

basedir=/project/daane/hussain/final_project/make_lastz_chains/obeta_ragtag_ncbi

csh $basedir/chain_obeta_cgob/cleanUp.csh

gunzip $basedir/chain_obeta_cgob/Cottoperca_gobio.Opsanus_beta.allfilled.chain.gz

/project/daane/shared/TOGA/toga.py $basedir/chain_obeta_cgob/Cottoperca_gobio.Opsanus_beta.allfilled.chain /project/daane/hussain/final_project/TOGA/cgob_ragtag.bed $basedir/chain_obeta_cgob/Cottoperca_gobio.2bit $basedir/chain_obeta_cgob/Opsanus_beta.2bit --pn toga_obeta_cgob -i /project/daane/hussain/final_project/TOGA/cgob_ragtag_isoforms.txt --nc /project/daane/shared/TOGA/nextflow_config_files --cesar_buckets 2,5,10,50
```

I got this error when I ran this script with cgob.bed and cgob_isoforms.txt files first:

```
cat cesar_jobs_crashed.txt
```


/project/daane/shared/TOGA/CESAR_wrapper.py ENSCGOT00000002507 38554 /project/daane/hussain/final_project/TOGA/obeta_ragtag_ncbi/toga_obeta_cgob/temp/toga_filt_ref_annot.hdf5 /project/daane/hussain/final_project/TOGA/obeta_ragtag_ncbi/toga_obeta_cgob/temp/genome_alignment.bst /project/daane/hussain/final_project/make_lastz_chains/obeta_ragtag_ncbi/chain_obeta_cgob/Cottoperca_gobio.2bit /project/daane/hussain/final_project/make_lastz_chains/obeta_ragtag_ncbi/chain_obeta_cgob/Opsanus_beta.2bit --cesar_binary /project/daane/shared/TOGA/CESAR2.0/cesar --uhq_flank 50 --temp_dir /project/daane/hussain/final_project/TOGA/obeta_ragtag_ncbi/toga_obeta_cgob/temp/cesar_temp_files --check_loss --alt_frame_del --memlim 2 CESAR JOB FAILURE       Input is corrupted! Reference sequence should start with ATG! Error! CESAR output is corrupted, target must start with ATG! Error! CESAR output is corrupted, target must start with ATG! Traceback (most recent call last):   File "/project/daane/shared/TOGA/CESAR_wrapper.py", line 2975, in <module>     realign_exons(cmd_args)   File "/project/daane/shared/TOGA/CESAR_wrapper.py", line 2940, in realign_exons     loss_report, del_mis_exons = inact_mut_check(                                  ^^^^^^^^^^^^^^^^   File "/project/daane/shared/TOGA/modules/inact_mut_check.py", line 1679, in inact_mut_check     split_stop_codons = detect_split_stops(                         ^^^^^^^^^^^^^^^^^^^   File "/project/daane/shared/TOGA/modules/inact_mut_check.py", line 1482, in detect_split_stops     position = exon_to_last_codon_of_exon[first_exon]                ~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^ KeyError: 1


So, I ran these commands so this transcript is not processed by TOGA:

```
sed '/ENSCGOT00000002507/d' cgob.bed > cgob_ragtag.bed
sed '/ENSCGOT00000002507/d' cgob_isoforms.txt > cgob_ragtag_isoforms.txt
```

## Processing TOGA Output Files

### Kidney Gene List

A list of kidney genes was made using Monarch Initiative website. The path to kidney genes is:
```
/project/daane/hussain/final_project/TOGA/kidney_genes.tsv
```

### Loss Summary

```
cd /project/daane/hussain/final_project/
mkdir toga_results
cd toga_results
mkdir loss_summ_data
cp /project/daane/hussain/final_project/TOGA/obeta_ncbi/toga_obeta_cgob/loss_summ_data.tsv obeta_ncbi.tsv
cp /project/daane/hussain/final_project/TOGA/obeta_ragtag_tama/toga_obet_cgob/loss_summ_data.tsv obeta_ragtag_tama.tsv
cp /project/daane/hussain/final_project/TOGA/obeta_ragtag_ncbi/toga_obeta_cgob/loss_summ_data.tsv obeta_ragtag_ncbi.tsv
touch toga_summary.py
python toga_summary.py obeta_ncbi.tsv,obeta_ragtag_ncbi.tsv,obeta_ragtag_tama.tsv
mv output.tsv obeta_summary.tsv
```
```
import sys
from collections import defaultdict
dataDict = defaultdict(dict)
specList=[]

# Check if command line arguments are provided
if len(sys.argv) < 2:
    print("Usage: python script.py file1.tsv,file2.tsv,file3.tsv")
    sys.exit(1)

# Get the comma-separated list of file paths from command line arguments
fileList = sys.argv[1].split(',')

for infile in fileList:
    spec=infile.split('.')[0]
    specList.append(spec)
    for line in open(infile):
        f=line.strip('\n').split('\t')
        if f[0]=="GENE":
            gene,result=f[1],f[2]
            dataDict[gene][spec]= result

with open('output.tsv', 'w') as outfile:
    header = "GENE\t" + "\t".join(specList)
    print(header, file=outfile)
    for gene in dataDict:
        row = gene
        for species in specList:
            result = dataDict[gene].get(species)
            if result is None:
                result = 'NA'
            row += "\t" + result
        print(row, file=outfile)
```

```
touch toga_filter_kidney_genes.py
python toga_filter_kidney_genes.py obeta_summary.tsv
mv test.txt obeta_kidney_summary.tsv
sed -i '/^$/d' obeta_kidney_summary.tsv
```
```
import sys
infile= sys.argv[1]
#Kidney Genes
geneList=['ENSCGOG00000000106','ENSCGOG00000000153',......./<too_long>/.......'ENSCGOG00000024645','ENSCGOG00000024660','ENSCGOG00000024687','ENSCGOG00000024740']
for line in open(infile):
    f=line.strip('\n').split('\t')
    if f[0]=='GENE':
        print(line, file=open('test.txt','a'))
    else:
        if f[0] in geneList:
            print(line,file=open('test.txt','a'))
```

```
touch count_common_status.py
python count_common_status.py obeta_summary.tsv
python count_common_status.py obeta_kidney_summary.tsv 
```
```
import sys

if len(sys.argv) < 2:
    print("Usage: python script.py <file_name>")
    sys.exit(1)

file_name = sys.argv[1]

same_in_2_and_3 = 0
same_in_2_and_4 = 0
same_in_3_and_4 = 0
same_in_all = 0

with open(file_name, 'r') as f:
    header = f.readline()
    for line in f:
        columns = line.strip().split('\t')
        col2, col3, col4 = columns[1], columns[2], columns[3]
        if col2 == col3:
            same_in_2_and_3 += 1
        if col2 == col4:
            same_in_2_and_4 += 1
        if col3 == col4:
            same_in_3_and_4 += 1
        if col2 == col3 == col4:
            same_in_all += 1

print(f"Number of genes with the same status in ncbi and ragtag_ncbi assemblies: {same_in_2_and_3}")
print(f"Number of genes with the same status in ncbi and ragtag_tama assemblies: {same_in_2_and_4}")
print(f"Number of genes with the same status in both ragtag assemblies: {same_in_3_and_4}")
print(f"Number of genes with the same status in all all assemblies: {same_in_all}")
```
```
touch count_values.py
python count_values.py obeta_summary.tsv
python count_values.py obeta_kidney_summary.tsv
```
```
import sys
from collections import defaultdict

file_name = sys.argv[1]

possible_values = ['I', 'L', 'UL', 'PI', 'M', 'PG', 'PM']

try:
    column_counts = defaultdict(lambda: defaultdict(int))

    # Open and read the file
    with open(file_name, 'r') as f:
        header = f.readline().strip().split('\t')
        num_columns = len(header)

        # Initialize counts for each value in each column
        for col_idx in range(num_columns):
            for value in possible_values:
                column_counts[col_idx][value] = 0

        # Count occurrences of each value in each column
        for line in f:
            columns = line.strip().split('\t')
            for col_idx, value in enumerate(columns):
                if value in possible_values:
                    column_counts[col_idx][value] += 1

    # Print count and percentage of each value in each column
    print("Count and percentage of each value in each column:")
    for col_idx, col_name in enumerate(header):
        total = sum(column_counts[col_idx].values())  # Total count for the column
        print(f"\n{col_name}:")
        for value in possible_values:
            count = column_counts[col_idx][value]
            percentage = (count / total * 100) if total > 0 else 0
            print(f"  {value}: {count} ({percentage:.2f}%)")
```

### Visualize Mutation

```
cd /project/daane/hussain/final_project/toga_results/
mkdir inact_mut_data
cd inact_mut_data/
cp /project/daane/hussain/final_project/TOGA/obeta_ncbi/toga_obeta_cgob/inact_mut_data.txt obeta_ncbi.txt
cp /project/daane/hussain/final_project/TOGA/obeta_ragtag_tama/toga_obet_cgob/inact_mut_data.txt obeta_ragtag_tama.txt
cp /project/daane/hussain/final_project/TOGA/obeta_ragtag_ncbi/toga_obeta_cgob/inact_mut_data.txt obeta_ragtag_ncbi.txt
```

Nephrin gene (_nphs1_)
```
/project/daane/shared/TOGA/supply/plot_mutations.py /project/daane/hussain/final_project/TOGA/cgob.bed obeta_ncbi.txt ENSCGOG00000009181 obeta_ncbi_nphs1.svg -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt
/project/daane/shared/TOGA/supply/plot_mutations.py /project/daane/hussain/final_project/TOGA/cgob.bed obeta_ragtag_tama.txt ENSCGOG00000009181 obeta_ragtag_tama_nphs1.svg -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt
/project/daane/shared/TOGA/supply/plot_mutations.py /project/daane/hussain/final_project/TOGA/cgob.bed obeta_ragtag_ncbi.txt ENSCGOG00000009181 obeta_ragtag_ncbi_nphs1.svg -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt
```

Podocin gene (_nphs2_)
```
/project/daane/shared/TOGA/supply/plot_mutations.py /project/daane/hussain/final_project/TOGA/cgob.bed obeta_ncbi.txt ENSCGOG00000011639 obeta_ncbi_nphs2.svg -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt
/project/daane/shared/TOGA/supply/plot_mutations.py /project/daane/hussain/final_project/TOGA/cgob.bed obeta_ragtag_tama.txt ENSCGOG00000011639 obeta_ragtag_tama_nphs2.svg -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt
/project/daane/shared/TOGA/supply/plot_mutations.py /project/daane/hussain/final_project/TOGA/cgob.bed obeta_ragtag_ncbi.txt ENSCGOG00000011639 obeta_ragtag_ncbi_nphs2.svg -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt
```

Gap Junction Alpha 1b gene (_gja1b_)
```
/project/daane/shared/TOGA/supply/plot_mutations.py /project/daane/hussain/final_project/TOGA/cgob.bed obeta_ncbi.txt ENSCGOG00000004743 obeta_ncbi_gja1b.svg -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt
/project/daane/shared/TOGA/supply/plot_mutations.py /project/daane/hussain/final_project/TOGA/cgob.bed obeta_ragtag_tama.txt ENSCGOG00000004743 obeta_ragtag_tama_gja1b.svg -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt
/project/daane/shared/TOGA/supply/plot_mutations.py /project/daane/hussain/final_project/TOGA/cgob.bed obeta_ragtag_ncbi.txt ENSCGOG00000004743 obeta_ragtag_ncbi_gja1b.svg -i /project/daane/hussain/final_project/TOGA/cgob_isoforms.txt
```

