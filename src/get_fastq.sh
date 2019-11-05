#! /bin/bash

# Hsapiens
data="/home/rstudio/disk"
mkdir -p $data
cd $data
mkdir -p sra_data
cd sra_data

#Liste des SRR d'intérêt
SRR="SRR3308952
SRR3308954
SRR3308956
SRR3308963
"

for srr in $SRR
do
# Télécharge les séquences correspondant aux SRR renseignés avant
fastq-dump $srr -O /home/rstudio/disk/sra_data -X 100000

done
