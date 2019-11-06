#! /bin/bash
## Run trimmomatic

cd /home/rstudio/disk

#Creation des dossiers pour stocker les fichiers nettoyés
mkdir data_trimmed
cd data_trimmed
mkdir paired
mkdir unpaired

cd ..

#Liste du numéro des patients
SRR="SRR3308952
SRR3308954
SRR3308956
SRR3308963
SRR3308971
SRR3308972
"

for fn in $SRR;
do
  #Trimmomatic sur des données Paired end
	java -jar /softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	  #On utilise que 7 coeurs sur les 8 de la machine
		-threads 7 \
		#Fichiers d'entrée
	  sra_data/$fn'_1.fastq' sra_data/$fn'_2.fastq' \
	  #Fichiers de sortie
	  data_trimmed/paired/$fn'_trimmed_paired_1.fastq' data_trimmed/unpaired/$fn'_trimmed_unpaired_1.fastq' \
	  data_trimmed/paired/$fn'_trimmed_paired_2.fastq' data_trimmed/unpaired/$fn'_trimmed_unpaired_2.fastq' \
	  #Fonctions à utiliser et arguments (les mêmes que dans le papier)
	  ILLUMINACLIP:adapateur.fa:2:30:10 \
	  LEADING:22 SLIDINGWINDOW:4:22 MINLEN:25

done

