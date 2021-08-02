# setwd("/ebio/abt6/yvoichek/smallproj/params/AraEnrich")
library(Matrix)
# AraEnrich will be a function that will find enrichment of genes in Arabidopsis vs a pre-defined group of genes.
# This infrastrcture will be built using Shiny - R library for building a gui
# we will start by coping the infrastructure we built in MATLAB:
# 1. We will do enrich vs gene-ontology
# 2. We will be able to plot the genes according to their location in the arabidopsis genome
# 3. We will enrich vs. all the predefined groups of genes 
# 4. We will be able to look at the list of genes and see their descriptions

#§ Loading the relevant libraries
# library(shiny)

AraEnrice_v1 <- function(genes, OG, bg_genes = NULL, n_res = 25) {
  genes2vec <- function(genes, n) {
    vec <- array(0, dim = c(1,n))
    vec[genes] <- 1  
    return(vec)
  }
  n_genes <- dim(OG$groups)[1];
  #############################################################################
  # If no background genes where specified, all genes will be considerd
  if(is.null(bg_genes)) {bg_genes <- 1:n_genes}
  bg_vec <- genes2vec(bg_genes, n_genes)
  #############################################################################
  
  #############################################################################
  # We want to caluclate a hypergeometric test on our genes vs. all the groups
  genes_vec <- genes2vec(intersect(genes, bg_genes), n_genes)
  # Size of the spaces - information for every group
  space_size <- bg_vec %*% OG$groups_bg
  space_size <- space_size[OG$g_info$background_index] 
  # size of sample genes
  sample_size <- genes_vec %*% OG$groups_bg
  sample_size <- sample_size[OG$g_info$background_index] 
  # size of groups
  groups_size <-  as.array(bg_vec %*% OG$groups)
  # size of intersect (interesting part)
  intersect_size <-  as.array(genes_vec %*% OG$groups)
  
  # Calculate the hypergeometric test
  res <- phyper(q=intersect_size-1,
                m = groups_size,
                n = space_size-groups_size,
                k = sample_size,
                log.p = TRUE, lower.tail = FALSE) / log(10)
  #############################################################################
  # Plot the results!
  sort_res <- sort(res, index.return = TRUE, decreasing = FALSE)
  output <- data.frame(name = I(array("",dim=c(1,n_res))),
                       log10_pval = I(array(NA,dim=c(1,n_res))),
                       group_index = I(array(NA,dim=c(1,n_res))),
                       groups_size = I(array(NA,dim=c(1,n_res))),
                       sample_size = I(array(NA,dim=c(1,n_res))),
                       intersect_size = I(array(NA,dim=c(1,n_res))),
                       set = I(array(NA,dim=c(1,n_res))))
  writeLines("using log10") 
  for(i in 1:n_res) {
    g_i <- sort_res$ix[i]
    output$group_index[i] <- g_i
    output$name[i] <- OG$g_info$name[g_i]
    output$log10_pval[i] <- round(sort_res$x[i],digits = 1)
    output$groups_size[i] <- groups_size[g_i]
    output$intersect_size[i] <- intersect_size[g_i]
    output$sample_size[i] <- sample_size[g_i]
    output$set[i] <- OG$g_info$set[g_i]
    
    rm(list = c('g_i'))
    
    writeLines(paste(
      i, "\t[", output$log10_pval[i] , "] [", 
      output$groups_size[i]," | ", output$sample_size[i], " | ", output$intersect_size[i], "]\t",
      output$set[i],':',output$name[i],'\r\n',
      sep=""))
  }
  invisible(output)
}


load("~/Desktop/forRobin/gene_infoV2.RData")
load("~/Desktop/forRobin/organize_groups_20191002.Rdata")


gene_names=cluster2_AAgene_colGO$V1
#gene_names=hightailGO$V5
genes <- c()
for(i in 1:length(gene_names)) {
  genes <- c(genes, which(gene_infoV2$Name == gene_names[i]))
}
#AraEnrice_v1(genes, OG,n_res=30)

background=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/Background.txt')

bg_gene_names <-  background$V1
bg_genes <- c()
for(i in 1:length(bg_gene_names)) {
  bg_genes <- c(bg_genes, which(gene_infoV2$Name == bg_gene_names[i]))
}
AraEnrice_v1(genes, OG,n_res=25,bg_genes)
##

#AraEnrice_v1(genes, OG,n_res=10, bg_genes)



############## For Paper for genes######


library(topGO)
library(KEGGREST)
library(FoldGO)
library(data.table)


background=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/Background.txt')

library("biomaRt")
#collect gene names from biomart
mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "athaliana_eg_gene",
                         host = 'plants.ensembl.org')
# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id",
                                        "go_id"), mart = mart)

#examine result
head (GTOGO)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
#make background my background
bg=GTOGO[GTOGO$ensembl_gene_id%in%background$V1,]


geneID2GO <- by(bg$go_id,
                bg$ensembl_gene_id,
                function(x) as.character(x))

geneNames <- names(geneID2GO)


aacluster1=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/Natural_AaSubgenomecluster1.txt',header=F)
atcluster1=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/Natural_AtSubgenomecluster1.txt',header=F)

aacluster2=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/Natural_AaSubgenomecluster2.txt',header=F)
atcluster2=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/Natural_AtSubgenomecluster2.txt',header=F)

atcluster3=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/Natural_AtSubgenomecluster3.txt',header=F)

lowtailtopgo=read.table(file='~/Desktop/lowtailGO.txt',header=F)
hightailtopgo=read.table(file='~/Desktop/hightailGO.txt',header=F)

ataafb=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/Asue_AT_Floralbud_AA_Floralbud_2FC.txt', header=F)  
ataar=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/Asue_AT_Rosette_AA_Rosette_2FC.txt', header=F)  

ataafb_atonly=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/AsueFloralBudLFC2_ATonly.txt', header=F)  
ataar_atonly=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/AsueRosetteLFC2_ATonly.txt', header=F)  

ataafb_aaonly=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/AsueFloralBudLFC2_AAonly.txt', header=F)  
ataar_aaonly=read.table(file='/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/htseq_1220/AsueRosetteLFC2_AAonly.txt', header=F)  

geneList <- factor(as.integer(geneNames %in% ataar$V1))
names(geneList) <- geneNames

GOdatar <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGor <- runTest(GOdatar, algorithm = "weight01", statistic = "fisher")
mysummaryr <- summary(attributes(resultTopGor)$score <= 0.001)
numsignifr <- as.integer(mysummaryr[[3]])
atr=GenTable(GOdatar, classic = resultTopGor,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignifr)

geneList <- factor(as.integer(geneNames %in% ataafb$V1))
names(geneList) <- geneNames

GOdatafb <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGofb <- runTest(GOdatafb, algorithm = "weight01", statistic = "fisher")
mysummaryfb <- summary(attributes(resultTopGofb)$score <= 0.001)
numsigniffb <- as.integer(mysummaryfb[[3]])
atfb=GenTable(GOdatafb, classic = resultTopGofb,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsigniffb)

geneList <- factor(as.integer(geneNames %in% ataar_atonly$V1))
names(geneList) <- geneNames

GOdatar_atonly <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGor_atonly <- runTest(GOdatar_atonly, algorithm = "weight01", statistic = "fisher")
mysummaryr_atonly <- summary(attributes(resultTopGor_atonly)$score <= 0.001)
numsignifr_atonly <- as.integer(mysummaryr_atonly[[3]])
atr_only=GenTable(GOdatar_atonly, classic = resultTopGor_atonly,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignifr_atonly)

geneList <- factor(as.integer(geneNames %in% ataafb_atonly$V1))
names(geneList) <- geneNames

GOdatafb_atonly <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGofb_atonly <- runTest(GOdatafb_atonly, algorithm = "weight01", statistic = "fisher")
mysummaryfb_atonly <- summary(attributes(resultTopGofb_atonly)$score <= 0.001)
numsigniffb_atonly <- as.integer(mysummaryfb_atonly[[3]])
atfb_only=GenTable(GOdatafb_atonly, classic = resultTopGofb_atonly,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsigniffb_atonly)

geneList <- factor(as.integer(geneNames %in% ataar_aaonly$V1))
names(geneList) <- geneNames

GOdatar_aaonly <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGor_aaonly <- runTest(GOdatar_aaonly, algorithm = "weight01", statistic = "fisher")
mysummaryr_aaonly <- summary(attributes(resultTopGor_aaonly)$score <= 0.001)
numsignifr_aaonly <- as.integer(mysummaryr_aaonly[[3]])
aaronly=GenTable(GOdatar_aaonly, classic = resultTopGor_aaonly,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignifr_aaonly)

geneList <- factor(as.integer(geneNames %in% ataafb_aaonly$V1))
names(geneList) <- geneNames

GOdatafb_aaonly <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGofb_aaonly <- runTest(GOdatafb_aaonly, algorithm = "weight01", statistic = "fisher")
mysummaryfb_aaonly <- summary(attributes(resultTopGofb_aaonly)$score <= 0.001)
numsigniffb_aaonly <- as.integer(mysummaryfb_aaonly[[3]])
aafbonly=GenTable(GOdatafb_aaonly, classic = resultTopGofb_aaonly,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsigniffb_aaonly)

R_plot=tableGrob(atr)
FB_plot=tableGrob(atfb)
atR_plot=tableGrob(atr_only)
aaR_plot=tableGrob(aaronly)
aaFB_plot=tableGrob(aafbonly)


library("cowplot")
tfirst_row=plot_grid(R_plot,FB_plot, labels = c('a', 'b'), align='l', ncol=2, nrow=1)
tsecond_row =plot_grid(atR_plot, aaFB_plot, labels = c('c', 'd'), align='l', ncol=2, nrow=1)
tthird_row=plot_grid(aaR_plot, ncol=1, nrow=1, labels='e')

pdf("~/Desktop/GO_tissues.pdf", width=20, height=20, useDingbats = F)
plot_grid(tfirst_row, tsecond_row, tthird_row, labels=c('', '',''), ncol=1,rel_heights = c(1.5, 0.3, 0.8))
dev.off()

pdf("~/Desktop/tableGO.pdf", width=12, height=20)
ggarrange(at1_plot,at2_plot,at3_plot,a1_plot,a2_plot,nrow=5,
          labels=c('a','b','c','d','e'))
dev.off()

geneList <- factor(as.integer(geneNames %in% lowtailtopgo$V1))
names(geneList) <- geneNames


GOdatalow <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGolow <- runTest(GOdatalow, algorithm = "weight01", statistic = "fisher")
mysummarylow <- summary(attributes(resultTopGolow)$score <= 0.001)
numsigniflow <- as.integer(mysummarylow[[3]])
atlow=GenTable(GOdatalow, classic = resultTopGolow,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsigniflow)

geneList <- factor(as.integer(geneNames %in% hightailtopgo$V1))
names(geneList) <- geneNames

GOdatahigh <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGohigh <- runTest(GOdatahigh, algorithm = "weight01", statistic = "fisher")
mysummaryhigh <- summary(attributes(resultTopGohigh)$score <= 0.001)
numsignifhigh <- as.integer(mysummaryhigh[[3]])
athigh=GenTable(GOdatahigh, classic = resultTopGohigh,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignifhigh)


 
geneList <- factor(as.integer(geneNames %in% atcluster1$V1))
names(geneList) <- geneNames

GOdata1 <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGo1 <- runTest(GOdata1, algorithm = "weight01", statistic = "fisher")
mysummary1 <- summary(attributes(resultTopGo1)$score <= 0.001)
numsignif1 <- as.integer(mysummary1[[3]])
at1_l=GenTable(GOdata1, classic = resultTopGo1,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif1)

geneList <- factor(as.integer(geneNames %in% atcluster2$V1))
names(geneList) <- geneNames

GOdata2 <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGo2 <- runTest(GOdata2, algorithm = "weight01", statistic = "fisher")
mysummary2 <- summary(attributes(resultTopGo2)$score <= 0.001)
numsignif2 <- as.integer(mysummary2[[3]])
at2_l=GenTable(GOdata2, classic = resultTopGo2,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif2)

geneList <- factor(as.integer(geneNames %in% atcluster3$V1))
names(geneList) <- geneNames

GOdata3 <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGo3 <- runTest(GOdata3, algorithm = "weight01", statistic = "fisher")
mysummary3 <- summary(attributes(resultTopGo3)$score <= 0.001)
numsignif3 <- as.integer(mysummary3[[3]])
at3_l=GenTable(GOdata3, classic = resultTopGo3,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif3)


geneList <- factor(as.integer(geneNames %in% aacluster1$V1))
names(geneList) <- geneNames

GOdata1a <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGo1a <- runTest(GOdata1a, algorithm = "weight01", statistic = "fisher")
mysummary1a <- summary(attributes(resultTopGo1a)$score <= 0.001)
numsignif1a <- as.integer(mysummary1a[[3]])
a1_l=GenTable(GOdata1a, classic = resultTopGo1a,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif1a)

geneList <- factor(as.integer(geneNames %in% aacluster2$V1))
names(geneList) <- geneNames

GOdata2a <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultTopGo2a <- runTest(GOdata2a, algorithm = "weight01", statistic = "fisher")
mysummary2a <- summary(attributes(resultTopGo2a)$score <= 0.001)
numsignif2a <- as.integer(mysummary2a[[3]])
a2_l=GenTable(GOdata2a, classic = resultTopGo2a,orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif2a)



# four dimension venn plot
pdf(file='~/Desktop/overlapcelldiv.pdf', height=2, width=4)
ggVennDiagram(x)
dev.off()

at1$color=rep("#66c2a5", nrow(at1))
at2$color=rep("#7570b3", nrow(at2))
at3$color=rep("#fc8d62", nrow(at3))
a1$color=rep("#e7298a", nrow(a1))
a2$color=rep("#e6ab02", nrow(a2))

at1$cluster=rep("At_cluster1", nrow(at1))
at2$cluster=rep("At_cluster2", nrow(at2))
at3$cluster=rep("At_cluster3", nrow(at3))
a1$cluster=rep("Aa_cluster1", nrow(a1))
a2$cluster=rep("Aa_cluster2", nrow(a2))


at1$fe=at1$Significant/at1$Expected
at2$fe=at2$Significant/at2$Expected
at3$fe=at3$Significant/at3$Expected
a1$fe=a1$Significant/a1$Expected
a2$fe=a2$Significant/a2$Expected

as_topGO=rbind(at1,at2,at3,a1,a2)

as_topGO$cluster=factor(as_topGO$cluster, levels=c('At_cluster1','At_cluster2','At_cluster3','Aa_cluster1','Aa_cluster2'))
as_topGO$fdr_pval=abs(-log10(as.numeric(as_topGO$classic)))


mybubble=ggplot(as_topGO, aes(y=fdr_pval, x=fe,color=color)) + 
  scale_size(range = c(.1, 15))+
  geom_point(aes(y=fdr_pval, x=fe, size=Significant,color=color,alpha=0.95)) +
  ylim(2,max(as_topGO$fdr_pval)+1) +
  xlim(-1,25) +
  geom_text(
    label=as_topGO$GO.ID, 
    nudge_x = 0.6, nudge_y = 0.6, 
    check_overlap = T
  ) +
  ylab("-log10(FDR)") +
  xlab("Fold Enrichment")  +
  scale_color_manual(name='Clusters',values=c('#66c2a5','#7570b3','#e7298a','#e6ab02','#fc8d62')) +
  theme(panel.background = element_blank(),plot.title = element_text(size=10,hjust = 0.5),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.title = element_text(size=10),  legend.text = element_text(size=10), axis.title=element_text(size=10)) +
  facet_grid(. ~ cluster, scales = "free", space = "free")


