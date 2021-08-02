rm(list=ls())
setwd('/Users/robin.burns/Google Drive/Suecica/RNAseq/SuecicaGenomeRNAmap/STAR')

knitr::opts_chunk$set(fig.width=3, fig.height=3, echo=FALSE, warning=FALSE, message=FALSE)
library(knitr)
library(limma)
library(edgeR)
library(locfit)
library(statmod)
library(Rcpp)
library(pheatmap)
library(gplots)
library(gdata)
library('RColorBrewer')
library(ape)
library(ggbiplot)
library(ggplot2)
library(reshape2)
library(ggExtra)
library(ggpubr)
library(dendextend)
library(vcd)
library("grid")
library(ggplotify)
library(pheatmap)
library(RColorBrewer)

tr_coord<-read.table('augustus_hints3sourcesAs_June2020.aug.renamed.gene.gff', header=F, stringsAsFactors=F)
tr_coord_AT<-tr_coord[which(tr_coord$V1 %in% c('Asue_scaffold1', 'Asue_scaffold2', 'Asue_scaffold3', 'Asue_scaffold4', 'Asue_scaffold5')),]

tr_coord_AA<-tr_coord[which(tr_coord$V1 %in% c('Asue_scaffold6', 'Asue_scaffold7', 'Asue_scaffold8', 'Asue_scaffold9',
                                               'Asue_scaffold10', 'Asue_scaffold11', 'Asue_scaffold12', 'Asue_scaffold13')),]


targets <- readTargets("mydescriptions.txt")
d <- readDGE(targets)
counts<-d$counts
counts=data.frame(counts)

#library sizes
barplot(d$samples$lib.size*1e-6, names=targets$description, ylab="Library size (millions)",cex.names=0.5,las=2)

#filter out genes with cross mapping Aar to Ath
countsAT<-counts[which(rownames(counts) %in% tr_coord_AT$V9),];
colnames(countsAT)<-c(targets$group)
testA<-countsAT[,which(colnames(countsAT)=='Aa')]
testT<-countsAT[,which(colnames(countsAT)=='At')]
testA=rowSums(testA)
testT=rowSums(testT)
pdf( "AthExp_logFCcutoff.pdf", width = 10, height = 7 )
hist(log(testT/testA),breaks=100)
tr=0
abline(v=tr, col='forest green')
dev.off()
qual<-c(); countsAT_filt<-c()
for (i in 1:dim(countsAT)[1]){
  if (testT[i]==0 & testA[i]>=quantile(testA)[4]){
    qual<-append(qual, 'bad')
  }else if (testT[i]==0 & testA[i]==0){
    qual<-append(qual, 'good')
    countsAT_filt<-rbind(countsAT_filt, countsAT[i,])
  }else if (log((testT)/(testA))[i]<0){
    qual<-append(qual, 'bad')
  }else{
    qual<-append(qual, 'good')
    countsAT_filt<-rbind(countsAT_filt, countsAT[i,])
  }
}
rownames(countsAT_filt)<-rownames(countsAT)[which(qual=='good')]


#filter out genes with cross mapping Ath to Aar

countsAA<-counts[which(rownames(counts) %in% tr_coord_AA$V9),];
colnames(countsAA)<-c(targets$group)
testA<-countsAA[,which(colnames(countsAA)=='Aa')]
testT<-countsAA[,which(colnames(countsAA)=='At')]
testT <- testT[,sample(1:45, size=20)]
testA=rowSums(testA)
testT=rowSums(testT)
pdf( "AarExp_logFCcutoff.pdf", width = 10, height = 7 )
hist(log(testA/testT),breaks=100)
tr=0
abline(v=tr, col='purple')
dev.off()
qual<-c(); countsAA_filt<-c()
for (i in 1:dim(countsAA)[1]){
  if (testA[i]==0 & testT[i]>=quantile(testT)[4]){
    qual<-append(qual, 'bad')
  }else if (testA[i]==0 & testT[i]==0){
    qual<-append(qual, 'good')
    countsAA_filt<-rbind(countsAA_filt, countsAA[i,])
  }else if (log((testA)/(testT))[i]<tr){
    qual<-append(qual, 'bad')
  }else{
    qual<-append(qual, 'good')
    countsAA_filt<-rbind(countsAA_filt, countsAA[i,])
  }
}
rownames(countsAA_filt)<-rownames(countsAA)[which(qual=='good')]

save.image(file = "Norm0920.RData")


##########################################################################################################
###homeologous pairs###

pairs<-read.table('onetoone_homeologs_withinAS.txt', header=T)

pairs<-pairs[,c('A.suecica.A.arenosa_sub.','A.suecica.A.thaliana_sub.')]
colnames(pairs)<-c('Asue_Aar','Asue_Ath')
pairs=pairs[!(duplicated(pairs$Asue_Ath)),]
pairs=pairs[!(duplicated(pairs$Asue_Aar)),]
pairs=pairs[!(duplicated(pairs)),]
#pairs=data.frame(pairs)

countsAT_filt1<-countsAT_filt[,-which(colnames(countsAT_filt)=='Aa')]
colnames(countsAT_filt1)<-targets$description[-which(colnames(countsAT_filt)=='Aa')]
countsAA_filt1<-countsAA_filt[,-which(colnames(countsAA_filt)=='At')]
colnames(countsAA_filt1)<-targets$description[-which(colnames(countsAA_filt)=='At')]

########reorganize data frames on both subgenomes to match rows by homeologs#####
library(dplyr)
pairs = pairs %>%
  arrange(Asue_Ath, Asue_Aar)

countsAT_filt1$Asue_Ath=rownames(countsAT_filt1)

countsAT_filt1 = countsAT_filt1 %>%
  arrange(Asue_Ath)

countsAA_filt1$Asue_Aar=rownames(countsAA_filt1)
countsAA_filt1 = countsAA_filt1 %>%
  arrange(Asue_Aar)


countsAT_filt_pairs=countsAT_filt1[countsAT_filt1$Asue_Ath%in%pairs$Asue_Ath,]
countsAA_filt_pairs=countsAA_filt1[countsAA_filt1$Asue_Aar%in%pairs$Asue_Aar,]


countsAT_filt_pairs=merge(countsAT_filt_pairs, pairs)
countsAT_filt_pairs=countsAT_filt_pairs[,c(1,95,2:94)]

countsAA_filt_pairs=merge(countsAA_filt_pairs, pairs)
countsAA_filt_pairs=countsAA_filt_pairs[,c(1,60,2:59)]


countsAT_filt_pairs = countsAT_filt_pairs %>%
  arrange(Asue_Ath, Asue_Aar)
countsAA_filt_pairs = countsAA_filt_pairs %>%
  arrange(Asue_Ath, Asue_Aar)

countsAT_filt_pairs_genes=countsAT_filt_pairs[,c(1,2)]
countsAA_filt_pairs_genes=countsAA_filt_pairs[,c(1,2)]
mymatchingpairs=merge(countsAT_filt_pairs_genes,countsAA_filt_pairs_genes, by=c('Asue_Ath','Asue_Aar'))


mymatchingpairs = mymatchingpairs %>%
  arrange(Asue_Ath, Asue_Aar)

countsAT_filt_pairs2=countsAT_filt_pairs[(countsAT_filt_pairs$Asue_Ath%in%mymatchingpairs$Asue_Ath & countsAT_filt_pairs$Asue_Aar%in%mymatchingpairs$Asue_Aar),]
countsAA_filt_pairs2=countsAA_filt_pairs[(countsAA_filt_pairs$Asue_Ath%in%mymatchingpairs$Asue_Ath & countsAA_filt_pairs$Asue_Aar%in%mymatchingpairs$Asue_Aar),]
countsAT_filt_pairs=countsAT_filt_pairs2
countsAA_filt_pairs=countsAA_filt_pairs2
save.image(file = "Norm0320.RData")

##############################################################################################################
#####normalize for gene length differences between Ath and Aar parts of the A. suecica genome##########

gff<-read.gff('augustus_Asuecica.gene.gff')
gff<-gff[which(gff$type=='gene'),]

ATgenes=countsAT_filt_pairs$Asue_Ath;AAgenes=countsAT_filt_pairs$Asue_Aar

normL<-c()
for (i in 1:length(ATgenes)){
  at<-ATgenes[i];aa<-AAgenes[i]
  Lat=abs(sum(gff[gff$attributes==at,'end']- gff[gff$attributes==at,'start']))
  Laa=abs(sum(gff[gff$attributes==aa,'end']- gff[gff$attributes==aa,'start']))
  normL<-rbind(normL, c(at,Lat,aa,Laa))
}
normL<-as.data.frame(normL)
colnames(normL)<-c('at','Lat','aa','Laa')
normL$mean = (as.numeric(as.character(normL$Lat)) + as.numeric(as.character(normL$Laa))) / 2

rownames(countsAT_filt_pairs)=countsAT_filt_pairs$Asue_Ath
rownames(countsAA_filt_pairs)=countsAT_filt_pairs$Asue_Aar
countsAT_filt_pairs=countsAT_filt_pairs[,c(3:ncol(countsAT_filt_pairs))]
countsAA_filt_pairs=countsAA_filt_pairs[,c(3:ncol(countsAA_filt_pairs))]

countsAT_filt_pairs_nl<-c();countsAA_filt_pairs_nl<-c();
for (i in 1:dim(countsAT_filt_pairs)[1]){
  countsAT_filt_pairs_nl<-rbind(countsAT_filt_pairs_nl, countsAT_filt_pairs[i,]*normL$mean[i]/as.numeric(as.character(normL$Lat[i])))
  countsAA_filt_pairs_nl<-rbind(countsAA_filt_pairs_nl, countsAA_filt_pairs[i,]*normL$mean[i]/as.numeric(as.character(normL$Laa[i])))
}

save.image(file = "Norm0320.RData")

##############################################################################################################
##################Normalize by effective library size in EdgeR##################

test<-paste(rep('At_', dim(countsAT_filt_pairs_nl)[2]), colnames(countsAT_filt_pairs_nl), sep='')
colnames(countsAT_filt_pairs_nl)<-test
groupAt<-gsub("_rep.","",test)
test<-paste(rep('Aa_', dim(countsAA_filt_pairs_nl)[2]), colnames(countsAA_filt_pairs_nl), sep='')
colnames(countsAA_filt_pairs_nl)<-test
groupAa<-gsub("_rep.","",test)

d<-DGEList(counts=cbind(countsAT_filt_pairs_nl, countsAA_filt_pairs_nl), group=c(groupAt, groupAa))


rownames(d) <- rownames(countsAT_filt_pairs_nl)

#########Filter out lowly expressed genes#############
keep <- d[((rowSums(cpm(d[,c(1:45)]) > 0) >=  10) & (rowSums(cpm(d[,c(46:69,72:82,85:87,92,93)]) > 0) >=  10) &
             (rowSums(cpm(d[,c(70,71,83,84)]) > 0) >=  1) &
             (rowSums(cpm(d[,c(88:91)]) > 0) >=  1) & (rowSums(cpm(d[,c(94:117,120:130,133:137,150,151)]) > 0) >=  10) & (rowSums(cpm(d[,c(118,119,131,132)]) > 0) >=  1) &
             (rowSums(cpm(d[,c(136:145)]) > 0) >=  2) & (rowSums(cpm(d[,c(146:148)]) > 0)>=  1)),]

keep <- keep[rowSums(cpm(keep)) > 5,]
df <- keep
df$samples$lib.size <- colSums(df$counts)


#######Normalize suecica library size by mean of two subgenomes########
mysuelibsize <- c()
for (i in 46:87) {
  suelibsize <- mean((sum(df$counts[,i]))+(sum(df$counts[,(i+48)])))
  mysuelibsize <- rbind(mysuelibsize, suelibsize)
}

df$samples$lib.size[46:87] <- mysuelibsize[1:42]
df$samples$lib.size[c(94:135)] <- mysuelibsize[1:42]

mysuelibsize <- c()
for (i in 88:93) {
  suelibsize <- mean((sum(df$counts[,i]))+(sum(df$counts[,(i+58)])))
  mysuelibsize <- rbind(mysuelibsize, suelibsize)
}


df$samples$lib.size[88:93] <- mysuelibsize[1:6]
df$samples$lib.size[c(146:151)] <- mysuelibsize[1:6]

###########Normalize by trimmed mean of means#################
df <- calcNormFactors(df, method='TMM')


############Get CPM and take mean of replicates###########
mygene_CPM <- cpm(df$counts, log=T, prior.count=1, normalized.lib.sizes=TRUE)
mygene_CPM <- data.frame(mygene_CPM)
colnames(mygene_CPM) = gsub("\\_rep1","",colnames(mygene_CPM))
colnames(mygene_CPM) = gsub("\\_rep2","",colnames(mygene_CPM))
colnames(mygene_CPM) = gsub("\\_rep3","",colnames(mygene_CPM))
colnames(mygene_CPM) = gsub("\\.","",colnames(mygene_CPM))
mygene_CPM.MEAN=t(apply(mygene_CPM, 1, function(x) tapply(x, colnames(mygene_CPM), mean)))

save.image(file = "Norm0320.RData")


##############################################################################################################
############Log fold change###########
###########homeologs##############
#Reorder
fr_aa=mygene_CPM.MEAN[,1:23]
fr_aa=data.frame(fr_aa)
fr_aa=fr_aa[,c(1,2,3,4,5:14,16,17,20,21,23,15,22,18,19)] #reorder 
fr_at=mygene_CPM.MEAN[,24:58]
fr_at=data.frame(fr_at)
fr_at=fr_at[,c(1:16,17:26,28,29,32,33,35,27,34,30,31)]


########Get LogF##########

t_as_aa <- fr_aa[,c(5:23)]
t_as_at <- fr_at[,c(17:35)]


t_aa_as <- fr_aa[,c(5:23)]
t_at_aa <- fr_at[,c(17:35)]

t_as_logFC=t_as_aa-t_as_at #subtracting logs is division


######Examine what is common in the tails of the logFC
t_as_logFC2=t_as_logFC[,c(1:15,18,19)] #ignore flower buds
quants <- c(0.05,0.95)
myquants=apply(t_as_logFC2 , 2 , quantile , probs = quants , na.rm = TRUE )

lowtail<-t_as_logFC2
for (i in 1:ncol(t_as_logFC2)) {
  lowtail[-which(t_as_logFC2[,i]<myquants[seq(1,34,2)][i]),i]<-0
}
lowtail=lowtail[apply(lowtail[,-1], 1, function(x) !all(x==0)),]
lowtail2<-t_as_logFC2
lowtail2=lowtail2[rownames(lowtail2)%in%rownames(lowtail),]
lowtail2=abs(lowtail2)
lowtailz=t(apply(lowtail2, 1, cal_z_score))


hightail<-t_as_logFC2
for (i in 1:ncol(t_as_logFC2)) {
hightail[-which(t_as_logFC2[,i]>myquants[seq(2,34,2)][i]),i]<-0
}
hightail=hightail[apply(hightail[,-1], 1, function(x) !all(x==0)),]
hightail2<-t_as_logFC2
hightail2=hightail2[rownames(hightail2)%in%rownames(hightail),]
hightailz=t(apply(hightail2, 1, cal_z_score))


###########################################################################
###############Species differences############################

#####homologs now, homeologs before#########


##############Species difference#########
#colors
#As=#762a83',Aa='#807dba',SynAa="#de77ae"
#As="#276419", At="#7fbc41",SynAt="#46afa9"


#Make PCA plot
library(ggbiplot)
#Arenosa subgenome
t_aa <- fr_aa[,c(1:19,22,23)]
speciesAA <-  c(rep("aa",4), rep("as", 15), rep("syn", 2))
speciesAA=factor(speciesAA, levels=c('aa', 'syn','as'))
PC_AA<-prcomp(t(t_aa));PCi<-data.frame(PC_AA$x,Species=speciesAA)
pc1=ggbiplot(PC_AA, group=speciesAA, ellipse = TRUE,circle = TRUE,var.axes = F) + theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.title = element_text(size=14),  legend.text = element_text(size=12), axis.title=element_text(size=14))+
  geom_point(aes(colour=speciesAA,size=3,alpha=0.8))+scale_color_manual(name="Species",values = c('#807dba',"#de77ae",'#762a83'), labels= c(expression(italic("A. arenosa"), "Synthetic",italic("A. suecica")))) + scale_size_continuous(name = "Density", guide = FALSE) + scale_alpha_continuous(name = "Density", guide = FALSE) + xlim(-2.5,2.5) + ylim(-2.5,4.5)

#Athaliana subgenome
t_at <- fr_at[,c(1:16,17:31,34,35)]
speciesAT <-  c(rep("at",16), rep("as", 15), rep("syn", 2))
speciesAT=factor(speciesAT, levels=c('at', 'syn','as'))
PC_AT<-prcomp(t(t_at));PCi<-data.frame(PC_AT$x,Species=speciesAT)
pc2=ggbiplot(PC_AT, group=speciesAT, ellipse = TRUE,circle = TRUE,var.axes = F) + theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.title = element_text(size=14),  legend.text = element_text(size=12), axis.title=element_text(size=14))+
  geom_point(aes(colour=speciesAT,size=3,alpha=0.8))+scale_color_manual(name="Species",values = c("#7fbc41","#46afa9","#276419"), labels= c(expression(italic("A. thaliana"),"Synthetic",italic("A. suecica")))) + scale_size_continuous(name = "Density", guide = FALSE) + scale_alpha_continuous(name = "Density", guide = FALSE) + xlim(-2.5,2.5) + ylim(-2.5,4.5)

######Significantly different genes expressed#######

#########Wilcoxon tests###########
#AA#
t_aa_pvalue=t_aa
t_aa_pvalue=as.matrix(t_aa_pvalue)
mypvalueaa=c()
library(stats)
for(i in 1:nrow(t_aa_pvalue)) {
  mygeneaa=t_aa_pvalue[i,]
  WT.t_aa_pvalue=wilcox.test(mygeneaa[c(1,2,3,4,20,21)],mygeneaa[c(5:19)]) #old columns before reps
  mypvalueaa_i=WT.t_aa_pvalue$p.value
  mypvalueaa=rbind(mypvalueaa,mypvalueaa_i)
}
mypvalueaa=as.data.frame(mypvalueaa)
rownames(mypvalueaa)=rownames(t_aa_pvalue)
mypvalueaa$genes=rownames(t_aa_pvalue);colnames(mypvalueaa)[1]="WilCoxProb"
mypvalueaa$fdr=p.adjust(mypvalueaa$WilCoxProb, method = "fdr", n = nrow(mypvalueaa))
#mypvalueaa_sig=mypvalueaa[mypvalueaa$WilCoxProb<0.05,]
#mypvalueaa$adj=p.adjust(mypvalueaa$WilCoxProb, method = "bonferroni", n = nrow(mypvalueaa))
mypvalueaa_sig=mypvalueaa[mypvalueaa$fdr<0.05,]

#AT#
t_at_pvalue=t_at
t_at_pvalue=as.matrix(t_at_pvalue)
mypvalueat=c()
for(i in 1:nrow(t_at_pvalue)) {
  mygeneat=t_at_pvalue[i,]
  WT.t_at_pvalue=wilcox.test(mygeneat[c(1:16,32,33)],mygeneat[c(17:31)])
  mypvalueat_i=WT.t_at_pvalue$p.value
  mypvalueat=rbind(mypvalueat,mypvalueat_i)
}
mypvalueat=as.data.frame(mypvalueat)
rownames(mypvalueat)=rownames(t_at_pvalue)
mypvalueat$genes=rownames(t_at_pvalue);colnames(mypvalueat)[1]="WilCoxProb"
mypvalueat$fdr=p.adjust(mypvalueat$WilCoxProb, method='fdr', n=nrow(mypvalueat))
mypvalueat_sig=mypvalueat[mypvalueat$WilCoxProb<=0.05,]
mypvalueat_sig=mypvalueat[mypvalueat$fdr<0.05,]

########Make Heatmap AA####
library(gplots)
library(RColorBrewer)

#######Aa
keep_taa3=t_aa[rownames(t_aa)%in%mypvalueaa_sig$genes,]  #####for text figure

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

keep_taa3 <- t(apply(keep_taa3, 1, cal_z_score))
library(dendextend)
AAmy_hclust_gene <- hclust(dist(keep_taa3), method = "ward.D2")
AAmy_gene_col <- cutree(tree = as.dendrogram(AAmy_hclust_gene), k = 3) #k=3
#three clusters
AAmy_gene_col=data.frame(cluster = ifelse(test = AAmy_gene_col == 1, (ifelse(test = AAmy_gene_col == 2,yes = "cluster2", no = "cluster3")), yes = "cluster1"))


AAmy_sample_col <- data.frame(species = rep(c("Aarenosa","Asuecica","Synthetic"), c(4,15,2)))
row.names(AAmy_sample_col) <- colnames(keep_taa3)

AAmy_colour = list(
  species = c(Aarenosa = "#807dba", Synthetic= "#de77ae",Asuecica = "#762a83"),
  cluster = c(cluster1 = "#e6ab02", cluster2 = "#e7298a", cluster3="#666666") # 
)

AAHeatmap=as.grob(pheatmap(keep_taa3,clustering_distance_rows='euclidean',clustering_distance_cols='euclidean',clustering_method="ward.D2",annotation_row = AAmy_gene_col, annotation_col = AAmy_sample_col,
         show_rownames=FALSE, annotation_colors = AAmy_colour,show_colnames=FALSE,cutree_cols=2,cutree_rows=3, color = colorRampPalette(rev(brewer.pal(n = 10,name ="RdBu")))(10)))


#####At 
keep_tat3=t_at[rownames(t_at)%in%mypvalueat_sig$genes,]  #from species test of signif

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

keep_tat3 <- t(apply(keep_tat3, 1, cal_z_score))


ATmy_hclust_gene <- hclust(dist(keep_tat3), method = "ward.D2")

ATmy_gene_col <- cutree(tree = as.dendrogram(ATmy_hclust_gene), k = 3)
ATmy_gene_col <- data.frame(cluster = ifelse(test = ATmy_gene_col == 1, ifelse(test = ATmy_gene_col == 2, yes='cluster2', no='cluster3'), yes='cluster1'))
                                                                        

ATmy_sample_col <- data.frame(species = rep(c("Athaliana","Asuecica","Synthetic"), c(16,15,2)))
row.names(ATmy_sample_col) <- colnames(keep_tat3)

ATmy_colour = list(
  species = c(Athaliana = "#7fbc41", Synthetic= "#46afa9",Asuecica = "#276419"),
  cluster = c(cluster1 = "#66c2a5", cluster2 = "#fc8d62", cluster3='#8da0cb')
)


ATHeatmap=as.grob(pheatmap(keep_tat3,clustering_distance_rows='euclidean',clustering_distance_cols='euclidean',clustering_method="ward.D2",annotation_row = ATmy_gene_col, annotation_col = ATmy_sample_col,
                           show_rownames=FALSE, annotation_colors = ATmy_colour,show_colnames=FALSE,cutree_cols=2,cutree_rows=3, color = colorRampPalette(rev(brewer.pal(n = 10,name ="RdBu")))(10)))


#############Get genes in the clusters for GO analysis###############

ATmy_gene_col2=ATmy_gene_col
ATmy_gene_col2$genes=rownames(ATmy_gene_col2)
cluster1_ATgene=ATmy_gene_col2[ATmy_gene_col2$cluster=='cluster1',]
cluster2_ATgene=ATmy_gene_col2[ATmy_gene_col2$cluster=='cluster2',]
cluster3_ATgene=ATmy_gene_col2[ATmy_gene_col2$cluster=='cluster3',]

cluster1_ATgene_colGO=myorthologs_Thit_filt_all[myorthologs_Thit_filt_all$V7%in%cluster1_ATgene$genes,]
write.table(cluster1_ATgene_colGO$V5, file='Natural_AtSubgenomecluster1.txt', row.names=F,col.names=F, quote=F, sep='\t')

cluster2_ATgene_colGO=myorthologs_Thit_filt_all[myorthologs_Thit_filt_all$V7%in%cluster2_ATgene$genes,]
write.table(cluster2_ATgene_colGO$V5, file='Natural_AtSubgenomecluster2.txt', row.names=F,col.names=F, quote=F, sep='\t')

cluster3_ATgene_colGO=myorthologs_Thit_filt_all[myorthologs_Thit_filt_all$V7%in%cluster3_ATgene$genes,]
write.table(cluster3_ATgene_colGO$V5, file='Natural_AtSubgenomecluster3.txt', row.names=F,col.names=F, quote=F, sep='\t')


AAmy_gene_col2=AAmy_gene_col
AAmy_gene_col2$genes=rownames(AAmy_gene_col2)
cluster3_AAgene=AAmy_gene_col2[AAmy_gene_col2$cluster=='cluster3',]
cluster1_AAgene=AAmy_gene_col2[AAmy_gene_col2$cluster=='cluster1',]
cluster2_AAgene=AAmy_gene_col2[AAmy_gene_col2$cluster=='cluster2',]

cluster1_AAgene_colGO=myorthologs_Thit_filt_all[myorthologs_Thit_filt_all$V7%in%cluster1_AAgene$genes,]
write.table(cluster1_AAgene_colGO$V5, file='Natural_AaSubgenomecluster1.txt', row.names=F,col.names=F, quote=F, sep='\t')

cluster2_AAgene_colGO=myorthologs_Thit_filt_all[myorthologs_Thit_filt_all$V7%in%cluster2_AAgene$genes,]
write.table(cluster2_AAgene_colGO$V5, file='Natural_AaSubgenomecluster2.txt', row.names=F,col.names=F, quote=F, sep='\t')

cluster3_AAgene_colGO=myorthologs_Thit_filt_all[myorthologs_Thit_filt_all$V7%in%cluster3_AAgene$genes,]
write.table(cluster3_AAgene_colGO$V5, file='Natural_AaSubgenomecluster3.txt', row.names=F,col.names=F, quote=F, sep='\t')


