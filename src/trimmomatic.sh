#! /bin/bash
## Run trimmomatic

cd /home/rstudio/disk

#Creation des dossiers pour stocker les fichiers nettoyés
mkdir -p data_trimmed
cd data_trimmed
mkdir -p paired
mkdir -p unpaired

cd ..

#Liste du numéro des patients
SRR="
SRR3308950
#SRR3308952
#SRR3308954
SRR3308958
SRR3308960
SRR3308961
SRR3308965
SRR3308967
SRR3308969
#SRR3308971
SRR3308977
SRR3308980
SRR3308984
#SRR3308956
#SRR3308963
#SRR3308972
SRR3308974
SRR3308978
SRR3308982
SRR3308976
"

for fn in $SRR;
do
  #Trimmomatic sur des données Paired end
  	  #On utilise que 7 coeurs sur les 8 de la machine
  	  #Fichiers d'entrée
  	  #Fichiers de sortie
  	  #Fonctions à utiliser et arguments (les mêmes que dans le papier)
  	  
	java -jar /softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
		-threads 7 \
	  sra_data/$fn'_1.fastq' sra_data/$fn'_2.fastq' \
	  data_trimmed/paired/$fn'_trimmed_paired_1.fastq' data_trimmed/unpaired/$fn'_trimmed_unpaired_1.fastq' \
	  data_trimmed/paired/$fn'_trimmed_paired_2.fastq' data_trimmed/unpaired/$fn'_trimmed_unpaired_2.fastq' \
	  ILLUMINACLIP:adaptateur.fa:2:30:10 \
	  LEADING:22 SLIDINGWINDOW:4:22 MINLEN:25

done

