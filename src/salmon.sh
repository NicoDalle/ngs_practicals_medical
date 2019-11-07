#! /bin/bash

cd /home/rstudio/disk
mkdir -p salmon

#Liste des SRR d'intérêt
SRR="SRR3308952
SRR3308954
SRR3308956
SRR3308963
SRR3308971
SRR3308972
"

## Run salmon

#Création de l'index à partir de la base de données, taille des k-mères : 25 (on a pas mal de reads assez petits)
#salmon index -t Hsap_cDNA.fa -i salmon/Hsap_index -k 25

for srr in $SRR :
do

#Quantification par salmon, -i: index, -l: type de librairie détecté automatiquement, -1 et -2: input
#-o: output, --validateMappings:
#--threads: nombre de coeurs dédiés
#--gcBias: option qui corrige s'il y a des changements dans le taux de GC
salmon quant -i salmon/Hsap_index -l A -1 data_trimmed/paired/$srr'_trimmed_paired_1.fastq' \
  -2 data_trimmed/paired/$srr'_trimmed_paired_2.fastq' \
  --validateMappings \
  -o salmon/$srr'_paired_quant' \
  --threads 7 --gcBias
echo $srr
done
