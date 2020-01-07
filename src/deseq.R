## DE analysis
library("tximport")
library("readr")
library(apeglm)
library("DESeq2",quietly = T)
library(ggplot2)

#On précise le dossier courant
setwd("/home/rstudio/disk/DESeq")

#Import des données, condition comporte les métadonnées, us comporte la table de compte produite par Salmon
condition <- read.table("condition.csv",header = T)
us=read.table("table_medical.csv")

# DE analysis : première comparaison du papier (avant/après chez les répondants)
#avec nos données, on prend en compte le facteur "temps" et "patient"
#Avec us[,condition$type=="responding"], on ne récupère que les niveaux d'expression des patients ayant répondu
dds <- DESeqDataSetFromMatrix(countData = round(us[,condition$type=="responding"],0),
                              colData = condition[condition$type=="responding",],
                              design = ~ patient + time)

#On enlève les gènes qui ne sont que très très peu exprimés, ceci permettra d'être moins stringent lors de l'analyse des gènes DE
#Pour les analyses, on met un FDR à 0,05 : 5% des gènes marqués comme DE ne le sont pas. Ceci est obtenu à partir de modèles statistiques
#En enlevant les gènes presque pas exprimés, on a un peu plus de gènes identifiés comme DE (mais le FDR est toujours le même)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#On définit la valeur "before" comme étant la référence, on comparera donc les conditions par rapport à la condition "before"
dds$time <- relevel(dds$time, ref = "before")
dds <- DESeq(dds)
dds_resp <- dds
resultsNames(dds)
#On shrink les données: on corrige le log change vis-à-vis de du fait que, quand le compte est faible, le log change est très variable comme il y a peu de comptes
#Ce qu'il faut saisir dans coef est une sortie de resultsNames(dds)
#(on fait juste ceci pour avoir le om qui correspond bien à ce qui est inscrit dans le fichier)
resLFC <- lfcShrink(dds, coef="time_after_vs_before", type="apeglm")
resLFC_resp <- resLFC

# 744 DE genes 5% FDR

#On fait le MA plot, pour vérifier la tête des données
plotMA(resLFC_resp,ylim=c(-4,4))
maplot <- ggplot(as.data.frame(resLFC_resp),aes(x=log10(baseMean),y=log2FoldChange,color=padj<0.05))+
  geom_point(mapping = aes(size=padj<0.05,alpha=padj<0.05,shape=padj<0.05,fill=padj<0.05))+theme_bw()+theme(legend.position = 'none')+
  scale_size_manual(values = c(0.1,1))+scale_alpha_manual(values = c(0.5,1))+
  scale_shape_manual(values=c(21,21))+scale_fill_manual(values = c("#999999","#05100e"))+
  scale_color_manual(values=c("#999999","#cc8167"))
maplot

#On ne sélectionne que les gènes qui sont exprimés différentiellement significativement
resLFC_resp$padj[is.na(resLFC_resp$padj)] <- 1
res_resp_sign <- resLFC_resp[resLFC_resp$padj<0.05,]
#Ici on classe les gènes dans l'ordre croissant (ou décroissant, avec le -) selon le fold change
#Ordre croissant, ce qui nous permet de voir les gènes les plus sous-exprimés
res_resp_sign <- res_resp_sign[order(res_resp_sign$log2FoldChange),]
gene_resp <- row.names(res_resp_sign)
#Ordre décroissant, ce qui nous permet de voir les gènes les plus sur-exprimés
res_resp_sign2 <- res_resp_sign[order(-res_resp_sign$log2FoldChange),]
gene_resp2 <- row.names(res_resp_sign2)

#PCA
vsd <- vst(dds_resp, blind=FALSE)
head(assay(vsd), 3)

#Plotter la PCA en basique, dans l'intgroup, on met tous les facteurs qu'on veut étudier, pour regarder s'ils ont un impact sur nos données
#Avec la table de métadonnées que l'on a, on peut mettre sex, time, type et patient
plotPCA(vsd, intgroup=c("sex", "time"))

#Plotter la PCA en joli
#Provient du tutoriel de DESeq, permet de mettre un facteur avec des couleurs, l'autre avec des formes
pcaData <- plotPCA(vsd, intgroup=c("sex", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type, shape=sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#On recommence les autres analyses du papier, de la même façon

# DE analysis : deuxième comparaison du papier (avant/après chez les non-répondants)
#avec nos données, on prend en compte le facteur "temps" et "patient"
dds <- DESeqDataSetFromMatrix(countData = round(us[,condition$type=="non_responding"],0),
                              colData = condition[condition$type=="non_responding",],
                              design = ~ patient + time)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$time <- relevel(dds$time, ref = "before") # Ici ce à quoi on va se comparer.
dds <- DESeq(dds)
dds_nonresp <- dds
res <- results(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="time_after_vs_before", type="apeglm")
resLFC_nonresp <- resLFC

# 33 DE genes 5% FDR
plotMA(resLFC_nonresp,ylim=c(-4,4))
maplot <- ggplot(as.data.frame(resLFC_nonresp),aes(x=log10(baseMean),y=log2FoldChange,color=padj<0.05))+
  geom_point(mapping = aes(size=padj<0.05,alpha=padj<0.05,shape=padj<0.05,fill=padj<0.05))+theme_bw()+theme(legend.position = 'none')+
  scale_size_manual(values = c(0.1,1))+scale_alpha_manual(values = c(0.5,1))+
  scale_shape_manual(values=c(21,21))+scale_fill_manual(values = c("#999999","#05100e"))+
  scale_color_manual(values=c("#999999","#cc8167"))
maplot
resLFC_nonresp$padj[is.na(resLFC_nonresp$padj)] <- 1
resLFC_nonresp[resLFC_nonresp$padj<0.05,]

#PCA

vsd <- vst(dds_nonresp, blind=FALSE)
head(assay(vsd), 3)

#Plotter la PCA en basique, dans l'intgroup, on met tous les facteurs qu'on veut étudier
plotPCA(vsd, intgroup=c("sex", "time"))

#Plotter la PCA en joli
pcaData <- plotPCA(vsd, intgroup=c("sex", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


# DE analysis : deuxième comparaison du papier (répondants/non répondants, avant traitement)
#avec nos données, on prend en compte le facteur "type" (i.e. est-ce que le patient a répondu ou on)
dds <- DESeqDataSetFromMatrix(countData = round(us[,condition$time=="before"],0),
                              colData = condition[condition$time=="before",],
                              design = ~ type)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$time <- relevel(dds$time, ref = "before") # Ici ce à quoi on va se comparer.
dds <- DESeq(dds)
dds_before <- dds
res <- results(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="type_responding_vs_non_responding", type="apeglm")
resLFC_bef <- resLFC

# 5 DE genes 5% FDR
plotMA(resLFC_bef,ylim=c(-4,4))
maplot <- ggplot(as.data.frame(resLFC_bef),aes(x=log10(baseMean),y=log2FoldChange,color=padj<0.05))+
  geom_point(mapping = aes(size=padj<0.05,alpha=padj<0.05,shape=padj<0.05,fill=padj<0.05))+theme_bw()+theme(legend.position = 'none')+
  scale_size_manual(values = c(0.1,1))+scale_alpha_manual(values = c(0.5,1))+
  scale_shape_manual(values=c(21,21))+scale_fill_manual(values = c("#999999","#05100e"))+
  scale_color_manual(values=c("#999999","#cc8167"))
maplot
resLFC_bef$padj[is.na(resLFC_bef$padj)] <- 1
resLFC_bef[resLFC_bef$padj<0.05,]

#PCA

#PCA

vsd <- vst(dds_before, blind=FALSE)
head(assay(vsd), 3)

#Plotter la PCA en basique, dans l'intgroup, on met tous les facteurs qu'on veut étudier
plotPCA(vsd, intgroup=c("sex", "type"))

#Plotter la PCA en joli
pcaData <- plotPCA(vsd, intgroup=c("sex", "time"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time, shape=sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

nous <- read.csv('nous.csv', sep = ' ', h =F)
intersect(nous, gene_resp)
 library(compare)
comparison <- compare(nous, gene_resp,allowAll=TRUE)
comparison$tM

a <- nous[nous %in% gene_resp]
