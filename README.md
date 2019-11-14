# NGS_Practicals_Medical

cd /home/rstudio/disk

Analyse des données de RNAseq des conditions 'Avant traitement chez les patients répondants' et 'Avant traitement chez les patients non répondants' du papier de Urup T. et al, 2017, BMC Cancer.

## Exécution des scripts présentés

# Exécution simple

Les scripts courts sont éxécutés directement dans le terminal en indicant le chemin menant vers les fichers .sh. Ex : pour éxécuter ./ngs_practicals_medical/src/get_fastq.sh, on entre simplement ./ngs_practicals_medical/src/get_fastq.sh dans le terminal.

# Exécution en fond (en nohup)

Les scripts peuvent prendre du temps pour être éxécutés. Les scripts sont lancés en fond via la fonction nohup, qui permet de pouvoir faire autre chose sur le terminal en parallèle et de ne pas afficher tout ce qui peut être retourné par les différentes fonctions sur le terminal. Ces informations sont cependant stockées dans un fichier nohup dans le répertoire courant.
Pour éxécuter ./ngs_practicals_medical/src/get_fastq.sh, il faut saisir 
nohup ./ngs_practicals_medical/src/get_fastq.sh [> <nom_du_fichier_de_stockage_des_infos_du_terminal>] &

## Téléchargement des données de séquençage de la banque de ARNm (plus ou moins)

Les données sont disponibles sur NCBI GEO au numéro GSE79671.
Les données pertinentes pour cette comparaison sont celles des tumeurs avant traitement :
  - non répondant 
"SRR3308950
SRR3308952
SRR3308954
SRR3308958
SRR3308960
SRR3308961
SRR3308965
SRR3308967
SRR3308969
SRR3308971
SRR3308977
SRR3308980
SRR3308984"
  - répondant 
"SRR3308956
SRR3308963
SRR3308972
SRR3308974
SRR3308978
SRR3308982
SRR3308976"

Pour chaque analyse/script, la liste des SRR est rentrée dans la valeur SRR. Cette liste sera parcourue par les boucles for à chaque fois, ce qui permet d'avoir les infos sur tous les patients.

Le téléchargement des données se fait depuis le script ./ngs_practicals_medical/src/get_fstq.sh. Les données brutes sont stockées sous la forme $srr_1.fastq (ou _2) dans le dossier sra_raw_data.

#Mise en forme des fichiers pour fastq et les quantifications

Dans le même script, les fichiers fastq sont modifiés pour que chaque ligne impaire (avec le nom des séquences) soit bien présentée pour pouvoir être lue par les programmes suivants.
Les données modifiées sont stockées comme $srr_1.fastq (ou_2) dans sra_data.

## Analyse de la qualité des reads avec fastqc

Pour analyser la qualité du séquençage, les fichiers sont analysés avec fastqc. Fastqc fait une analyse globale de la qualité des reads, de la présence de séquences surreprésentées... Le script correspondant est ./ngs_practicals_medical/src/fastqc.sh, et est appliqué aux données présentes dans sra_data. On obtient un rapport par fichier SRR. Un rapport global est établi avec multiqc (script : ./ngs_practicals_medical/src/multiqc.sh). Toutes les données de fastqc et multiqc sont stockées dans data_qc.

Rq : le script fastqc_trimmed.sh permet de faire exactement la même chose, mais il est paramétré pour choisir les données nettoyées par trimmomatic (cf plus bas).

## Nettoyage des données

Le nettoyage des données de séquençage est réalisé par trimmomatic (script ./ngs_practicals_medical/src/trimmomatic.sh). Pour chaque read, trimmomatic enlève les bases dont la confiance est très faible (avec la fonction SLIDINGWINDOW -> enlève les bases d'une fen^tre glissante si la confiance moyenne de ces bases est en dessous d'une certaine valeur) et les séquences qui ont été trouvées anormalement surreprésentées par fastqc (avec la fonction ILLUMINACLIP). 
Dans notre cas, seulement les séquences correspondant aux séquences polyA/polyT et les aaptateurs Illumina ont été ajoutés dans le fichier adaptateur.fa.
Les reads nettoyés qui sont en dessous d'une certaine taille sont enelvés avec la fonction MINLEN et les bases du début du read dont la qualité est inféreure à un certain seuil sont enlevées avec LEADING (prises une par une cette fois, contrairement à SLIDINGWINDOW). 
Cependant, il se peut que seulement une des deux séquences correpondant au même read ne soit de qualité suffisante pour être conservée. Dans ce cas, ces séquences nettoyées (qui ne sont plus pairées) sont renvoyées dans les fichiers ./data_trimmed/unpaired/$srr_trimmed_unpaired_1.fatsq (ou_2). Ces séquences sont orphelines et ne seront pas utilisées pour la suite. Les séquences nettoyées mais ayant toujours leur séquence correspondante à l'autre bout du read sont renvoyées dans ./data_trimmed/paired/$srr_trimmed_paired_1.fatsq (ou_2).

Les paramètres utilisés ici sont les mêmes que ceux de l'article (sauf la fonction HEADCROP qui a été enlevée).

# Analyse du nettoyage

Le contrôle qualité des séquences nettoyées est réalisé avec fastqc et multiqc (scripts ./ngs_practicals_medical/src/fastqc_trimmed.sh et ./ngs_practicals_medical/src/multiqc_trimmed.sh) comme décrit précédemment.

## Quantification des transcrits et mapping sur le transcriptome humain de référence

# Téléchargement du transcriptome humain

Le transcriptome humain est récupéré depuis ensembl.org/biomart. La version GRCh38.p13 de Ensembl Gene 98 a été utilisée. Le script utilisé est ./ngs_practicals/src/get_cdna.sh. Le transcriptome est enregistré en tant que Hsap_cDNA.fa

# Quantification et mapping des reads sur le trnascriptome avec salmon

Pour cette étape, la fonction salmon a été utilisée à partir des données nettoyées et du transcriptome téléchargé juste avant. Le script est ./ngs_practicals_medical/src/salmon.sh.
Rq : salmon permet de quantifier des reads de RNAseq sur une banque de cDNA seulement. pour une quantification sur un génome entier, il faut utiliser STAR par exemple (cf plus bas).
Salmon fonctionne en 2 étapes. D'abord la création d'un index à partir du transcriptome dans le dossier salmon/Hsap_index. Cet index n'a pas besoin d'être regénéré une fois fait. 
Ensuite, salmon prend les 2 fichiers de séquences de séquencçage nettoyées et appariées et les mappe sur l'index et quantifie le nombre de reads qui se réfèrent au même transcrit. Le type de librairie n'étant pas précisé dans le papier, la foncion -l A est mise, ce qui permet à salmon de voir quel type de librairie est le plus probable d'avoir été utilisé par les auteurs.
Pour chaque paire de fichiers de séquences, les données sont stockées dans salmon/$srr_paired_quant. 

# Résultats

Le fichier cmd_info.json fait un récapitulatif de l'analyse qui a été faite. La table est dans quant.sf. Ce fichier comprend plusieurs colonnes : 
  - Name : référence du gène, référence du transcrit, nom du gène
  - Length : longueur du transcrit
  - EffectiveLength : ?
  - TPM : transcrits par million : proportion du transcrit parmi les ARN analysés. Ceci permet de prendre en compte le niveau d'expression du gène, la longueur du gène (plus il est long, plus on a de chance d'avoir de reads dessus) et la profondeur du séquençage. Ceci reflète le niveau d'expression normalisé par rapport à la longueur du transcrit et de la profondeur du séquençage. Pas forcément pertinent pour notre étude : ne varie pas si on a une augmentation globale de la quantité d'ARN et le TPM d'un transcrit dont l'expression ne change pas peut changer si tous les autres transcrits sont différemment exprimés.
  - NumReads : Nombre de reads mappés sur ce transcrit. C'est ce que l'on va utiliser par la suite (il faudra cependant corriger ceci par la profodeur du séquençage). Pas besoin de normalisation par la taille du transcrit comme on va comparer l'expression d'un même transcrit dans 2 conditions différentes.
Enfin, durant que salmon s'exécute, il renvoie les informations de la qualité du mapping en pourcentage pour chaque paire de fichiers de séquences. Attention à bien conserver le nohup quand salmon.sh est exécuté !! Cette qualité du mapping est en fait le nombre de reads que salmon a réussi à mapper sur le transcriptome de référence. Dans notre cas, la qualité du mapping est très faible (entre 15 et 30%), bien en dessous d'une valeur acceptable (au moins 80% dans notre cas : RNAseq sur un transcriptome d'une espèce dont le génome est bien connu et annoté). 

## Quantification et mapping des reads sur le génome humain