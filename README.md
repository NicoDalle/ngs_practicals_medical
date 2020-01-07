# NGS_Practicals_Medical

cd /home/rstudio/disk

Analyse des données de RNAseq des conditions 'Avant traitement chez les patient.e.s répondant.e.s' et 'Avant traitement chez les patient.e.s non répondant.e.s' du papier de Urup T. et al, 2017, BMC Cancer.

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

Pour chaque analyse/script, la liste des SRR est rentrée dans la valeur SRR. Cette liste sera parcourue par les boucles for à chaque fois, ce qui permet d'avoir les infos sur tous les patient.e.s.

Le téléchargement des données se fait depuis le script ./ngs_practicals_medical/src/get_fstq.sh. Les données brutes sont stockées sous la forme $srr_1.fastq (ou _2) dans le dossier sra_raw_data.

#Mise en forme des fichiers pour fastq et les quantifications

Dans le même script, les fichiers fastq sont modifiés pour que chaque ligne impaire (avec le nom des séquences) soit bien présentée pour pouvoir être lue par les programmes suivants.
Les données modifiées sont stockées comme $srr_1.fastq (ou_2) dans sra_data.

## Analyse de la qualité des reads avec fastqc

Pour analyser la qualité du séquençage, les fichiers sont analysés avec fastqc. Fastqc fait une analyse globale de la qualité des reads, de la présence de séquences surreprésentées... Le script correspondant est ./ngs_practicals_medical/src/fastqc.sh, et est appliqué aux données présentes dans sra_data. On obtient un rapport par fichier SRR. Un rapport global est établi avec multiqc (script : ./ngs_practicals_medical/src/multiqc.sh). Toutes les données de fastqc et multiqc sont stockées dans data_qc.

Rq : le script fastqc_trimmed.sh permet de faire exactement la même chose, mais il est paramétré pour choisir les données nettoyées par trimmomatic (cf plus bas).

## Nettoyage des données

Le nettoyage des données de séquençage est réalisé par trimmomatic (script ./ngs_practicals_medical/src/trimmomatic.sh). Pour chaque read, trimmomatic enlève les bases dont la confiance est très faible (avec la fonction SLIDINGWINDOW -> enlève les bases d'une fenêtre glissante si la confiance moyenne de ces bases est en dessous d'une certaine valeur) et les séquences qui ont été trouvées anormalement surreprésentées par fastqc (avec la fonction ILLUMINACLIP). 
Dans notre cas, seulement les séquences correspondant aux séquences polyA/polyT et les adaptateurs Illumina ont été ajoutés dans le fichier adaptateur.fa.
Les reads nettoyés qui sont en dessous d'une certaine taille sont enlevés avec la fonction MINLEN et les bases du début du read dont la qualité est inféreure à un certain seuil sont enlevées avec LEADING (prises une par une cette fois, contrairement à SLIDINGWINDOW). 
Cependant, il se peut que seulement une des deux séquences correpondant au même read ne soit de qualité suffisante pour être conservée. Dans ce cas, ces séquences nettoyées (qui ne sont plus pairées) sont renvoyées dans les fichiers ./data_trimmed/unpaired/$srr_trimmed_unpaired_1.fatsq (ou_2). Ces séquences sont orphelines et ne seront pas utilisées pour la suite. Les séquences nettoyées mais ayant toujours leur séquence correspondante à l'autre bout du read sont renvoyées dans ./data_trimmed/paired/$srr_trimmed_paired_1.fatsq (ou_2).

Les paramètres utilisés ici sont les mêmes que ceux de l'article (sauf la fonction HEADCROP qui a été enlevée).

# Analyse du nettoyage

Le contrôle qualité des séquences nettoyées est réalisé avec fastqc et multiqc (scripts ./ngs_practicals_medical/src/fastqc_trimmed.sh et ./ngs_practicals_medical/src/multiqc_trimmed.sh) comme décrit précédemment.

## Quantification des transcrits et mapping sur le transcriptome humain de référence

# Téléchargement du transcriptome humain

Le transcriptome humain est récupéré depuis ensembl.org/biomart. La version GRCh38.p13 de Ensembl Gene 98 a été utilisée. Le script utilisé est ./ngs_practicals/src/get_cdna.sh. Le transcriptome est enregistré en tant que Hsap_cDNA.fa

# Quantification et mapping des reads sur le transcriptome avec salmon

Pour cette étape, la fonction salmon a été utilisée à partir des données nettoyées et du transcriptome téléchargé juste avant. Le script est ./ngs_practicals_medical/src/salmon.sh.
Rq : salmon permet de quantifier des reads de RNAseq sur une banque de cDNA seulement. Pour une quantification sur un génome entier, il faut utiliser STAR par exemple (cf plus bas).
Salmon fonctionne en 2 étapes. D'abord la création d'un index à partir du transcriptome dans le dossier salmon/Hsap_index. Cet index n'a pas besoin d'être regénéré une fois fait. 
Ensuite, salmon prend les 2 fichiers de séquences de séquençage nettoyées et appariées et les aligne sur l'index et quantifie le nombre de reads qui se réfèrent au même transcrit. Le type de librairie n'étant pas précisé dans le papier, la foncion -l A est mise, ce qui permet à salmon de voir quel type de librairie est le plus probable d'avoir été utilisé par les auteurs.
Pour chaque paire de fichiers de séquences, les données sont stockées dans salmon/$srr_paired_quant. 
/!\ Bien garder le fichier nohup, car il contient le pourcentage de mapping des reads sur le transcriptome de référence, qui normalement doit être au dessus de 70%... Dans notre cas, on tournait plus autour de 15-30%...

# Résultats

Le fichier cmd_info.json fait un récapitulatif de l'analyse qui a été faite. La table est dans quant.sf. Ce fichier comprend plusieurs colonnes : 
  - Name : référence du gène, référence du transcrit, nom du gène
  - Length : longueur du transcrit
  - EffectiveLength : la longueur effective du transcrit telle que modélisée par salmon, à partir de la dsitrbution de la taille des fragments, des biais liés à la séquence et à la teneur en gc. (Ne sera pas utile pour nous...)
  - TPM : transcrits par million : proportion du transcrit parmi les ARN analysés. Ceci permet de prendre en compte le niveau d'expression du gène, la longueur du gène (plus il est long, plus on a de chance d'avoir de reads dessus) et la profondeur du séquençage. Ceci reflète le niveau d'expression normalisé par rapport à la longueur du transcrit et de la profondeur du séquençage. Pas forcément pertinent pour notre étude : ne varie pas si on a une augmentation globale de la quantité d'ARN et le TPM d'un transcrit dont l'expression ne change pas peut changer si tous les autres transcrits sont différemment exprimés.
  - NumReads : Nombre de reads mappés sur ce transcrit. C'est ce que l'on va utiliser par la suite (il faudra cependant corriger ceci par la profodeur du séquençage). Pas besoin de normalisation par la taille du transcrit comme on va comparer l'expression d'un même transcrit dans 2 conditions différentes.
Enfin, pendant que salmon s'exécute, il renvoie les informations de la qualité du mapping en pourcentage pour chaque paire de fichiers de séquences. Attention à bien conserver le nohup quand salmon.sh est exécuté !! Cette qualité du mapping est en fait le nombre de reads que salmon a réussi à mapper sur le transcriptome de référence. Dans notre cas, la qualité du mapping est très faible (entre 15 et 30%), bien en dessous d'une valeur acceptable (au moins 80% dans notre cas : RNAseq sur un transcriptome d'une espèce dont le génome est bien connu et annoté). 

Pour essayer de comprendre ce qui a pu se passer lors de la préparation de la librairie, nous avons décidé d'aligner les reads sur le génome humain. Ceci devrait nous permettre d'y voir plus clair. Ceci ne peut pas être fait à partir de Salmon, mais nous allons utiliser un programme qui fonctionne de manière assez similaire sur le principe : STAR.

## Quantification et mapping des reads sur le génome humain

#Téléchargement du génome de référence et de son annotation

La 38e version du génome humain est téléchargée depuis le site de l'UCSC par le scirpt get_genome.sh. En parallèle, l'annotation de ce génome est aussi téléchargée, ce qui permet d'avoir toutes les informations concernant les séquences sur lesquelles vont se mapper nos reads. Ceci est fait avec le même script.

#Quantification et mapping des reads sur le génome avec STAR

Le script permettat de faire fonctionner STAR s'appelle star.sh. Il utilise les fichiers .fastq nettoyés après trimmomatic et le génome et son annotation. Dans un premier temps, un index du génome humain est généré avec --runMode genomeGenerate à partir du génome et de son annotation. Malheureusement, nous n'avons pas pu générer l'index par nous-mêmes (même en ne prenant en compte que les chromosomes entiers, ces "améliorations" sont présentées commentées dans le code star.sh). L'index généré par une autre machine a été utilisé, se trouvant dans ./Corentin/hg38_genome. Le mapping et la quantification ont ensuite été réalisées en même temps par STAR en utilisant nos fichiers nettoyés et l'index. Les fichiers de sortie sont sous format .sam, ils sont compressés en .bam avec Samtools.

#Vérification du mapping

Les analyses du mapping sur le génome humain dans sa globalité sont faites pour chaque échantillon avec qualimap. On se sert des fichiers .bam renvoyés par STAR, qu'il faut d'abord classer par ordre alphabétique de gène, et pas par ordre de position sur le chromosome (comme ils sont classés par STAR initialement). Une fois ceci fais, on peut appliquer qualimap, ce qui dit nous donner des informations sur la qualité du mapping des reads sur le génome de référence qui a été fait par STAR. En théorie, on devrait avoir un très bon mapping. 
Cependant, nous n'avons pas pu obtenir de résultats de qualimap suite à des problèmes de versios de JAVA qui n'étaient pas à jour... 

#Résultats du mapping

Les résultats du mapping par STAR auraient pu être analysés si nous avions eu qualimap... Du coup nous ne savons pas exactement ce qui a pu se passer, mais nous pensons que la dépletion en ARNr ne s'est pas bien faite, ce qui a conduit à une grosse contamination par les ARNr.

##Analyse statistique de la table de compte obtenue par Salmon

Toute l'aalyse statistique est présentée dans un seul script R : deseq.R. Il y a toutes les étapes de l'analyse pour chacune des comparaisons.

#Récupération des métadonnées de chaque patient.e

Pour chaque échantillon, le NCBI donne des informations sur le sexe du patient.e, sa réponse, le numéro de lot de séquençage duquel il fait partie... Ce sont les métadonnées qui sont récupérées directement sur GEO et stockées dans le fichier condition.csv.

Pour les analyses statistiques, nous sommes reparti.e.s chacun.e.s avec l'ensemble des quantifications pour tous les patient.e.s, ce qui nous a permis de faire toutes les comparaisons du papier.

#Analyse d'expression différentielle avec DESeq

QUELS PARAMETRES ONT ETE PRIS EN COMPTE

Pour chaque comparaison, les données des patient.e.s pertinents sont récupérées et rentrées dans DESeq. Les gènes présentant très peu de reads (moins de 10 sur la totalité des patient.e.s, ce chiffre est arbitraire mais permet juste d'éliminer les gènes qui ne sont pas détectés). Ceci nous permet de gagner un peu de puissance lors de l'analyse statistique par DESeq. Ce test nous permet de savoir si le gène est différentiellement exprimé entre le témoin et la condition à analyser. On met un FDR (False Discovery Rate) à 0,05, ce qui veut dire que, d'après des modèles statistiques, 5 % des gènes marqués comme différentiellement exprimés ne le sont pas.

#MA Plot

On vérifie ensuite que les données ont une distribution attendue en traçant le MA plot des données brute set shrinkées (càd corrigées pour les petits comptes, ce qui a tendance à faire exploser la variance). Le MA plot est un graphe qui présente le logarithme de ratio entre les deux conditions en fonction du nombre de compte pour chaque gène. On peut détecter les gènes différentiellement exprimés.

#Récupération et identification des gènes différentiellement exprimés

On récupère les gènes dont la p-value est iférieure à 0,05. On peut ensuite récupérer les gènes sur- et sous-exprimés en triant simplement la liste. Les 5 gènes les plus sur- et sous-exprimés ont été idetifiés en utilisant les identifiants ENSEMBL manuellement. Je pense qu'il doit être possible de pouvoir faire cette identification (ainsi qu'une analyse de Gene Ontology) informatiquement mais nous n'avons pas eu le temps de la réaliser.

#PCA

Pour vérifier si les autres facteurs autre que celui qui est intéressant pour la comparaison jouent un rôle sur la structure des données, des Analyses par Composantes Principales ont été faites en utilisant différentes combinaisons de paramètres (sexe, type, time, patient). Les échantillons venant du même patient.e sont proches. De plus, on peut remarquer des clusters qui ne sont pas expliqués par les facteurs autres. Cependant, le facteur qui nous intéressait n'était pas vraiment responsable du clustering des données (en tout cas au moins pour la comparaison 3 : répondant.e.s et non-répondat.e.s avant triatement, ceci paraît cohérent avec le fait que seulement 5 gènes sont différentiellement exprimés entre ces deux conditions.)

Les analyses sont répétées piur chacune des comparaisons.

##Suite du travail

A partir de la liste des gènes qui sont différentiellement exprimés, il serait intéressant de les identifier et de vois si certains processus biologiques sont touchés par le traitement. Une analyse de Gene Ontology serait utile dans ce cas. Le papier a aussi essayé de reformer des réseaux de gènes interagissant entre eux à partir des données de séquençage. Il pourrait être intéressant de vois si le réseau obtenu dans le papier colle à ce que nous pourrions avoir avec notre analyse.