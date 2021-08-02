#### HMM for genotyping of samples for Asuecica Genetic Map####
#From Daniele#

#Make genetic map filtering for MQ 60 and appropriate number of reads per snp in mpileup
#min 2 for F2s
setwd("/Users/robin.burns/Documents/002Asuecica/016GeneticMap")

library(data.table)
library(genoPlotR)
library(qtl)
library(ggplot2)

#######################Change this section accoring to chromosome number used here Chr8###################
P1<-read.table('ASS3a.mpileup.8', header=F, na.strings="", fill=T)[,1:5]
P2<-read.table('AS150.mpileup.8', header=F, na.strings="", fill=T)[,1:5]
parents<-merge(P1,P2, by="V2")[,c(2,1,5,9)]
colnames(parents)<-c('CHROM', 'POS','ASS3a', 'AS150')
parents_new<-c()
for (i in 1:dim(parents)[1]){
  my<-table(unlist(strsplit(toupper(parents$ASS3a[i]), split='')))
  gen_p1<-names(my)[which(my==max(my))]
  my<-table(unlist(strsplit(toupper(parents$AS150[i]), split='')))
  gen_p2<-names(my)[which(my==max(my))]
  if (length(gen_p1)==1 & length(gen_p2)==1){
    parents_new<-rbind(parents_new, cbind(as.vector(parents[i,c(1,2)]),cbind(gen_p1, gen_p2)))
  }
}
parents_new<-parents_new[as.character(parents_new$gen_p1)!=as.character(parents_new$gen_p2),]
parents_new<-parents_new[as.character(parents_new$gen_p2)!="*" & as.character(parents_new$gen_p2)!="N",]

parents_new$gen_p1=gsub("\\,", "\\.", parents_new$gen_p1, perl = T)
parents_new$gen_p2=gsub("\\,", "\\.", parents_new$gen_p2, perl = T)
parents_new<-parents_new[as.character(parents_new$gen_p1)!=as.character(parents_new$gen_p2),]
parents_new<-parents_new[as.character(parents_new$gen_p1)==".",]
#
#parents_new=parents_new[,c(1,2,3,5)];colnames(parents_new)[3]="gen_p1"

chr_files<-list.files('/Users/robin.burns/Documents/002Asuecica/016GeneticMap', pattern='[.]8$')
chr_files<-chr_files[grep('Merged', chr_files)]
for (file in chr_files){
  my<-read.table(paste('./', file, sep=''), na.strings="", header=F, fill=T)[,c(2,5)]
  my<-my[which(my$V2 %in% parents_new$POS),]
  colnames(my)<-c('POS','gen')
  my$gen=gsub("\\,", "\\.", my$gen, perl=T)
  test<-merge(parents_new, my, by='POS', all=T)
  gen<-toupper(test$gen)
  parents_new$gen<-as.vector(gen)
  colnames(parents_new)<-append(colnames(parents_new)[1:length(colnames(parents_new))-1],gsub("\\.mpileup.8", "", file, perl = T))
}

parents_good<-parents_new
for (i in 1:dim(parents_new)[1]){
  my<-parents_new[i,5:ncol(parents_new)]
  gen<-c()
  for (j in 1:length(my)){
    if (is.na(my[j])){
      gen<-append(gen,NA)
    }else{
      p1<-length(which(unlist(strsplit(as.character(my[j]), split=''))==as.character(parents_new$gen_p1[i]))) #this is wrong as it was j fixed to i
      p2<-length(which(unlist(strsplit(as.character(my[j]), split=''))==as.character(parents_new$gen_p2[i])))
      if (p1!=0 & p2!=0){
        gen<-append(gen,'0.5')
      } else{
        if (p1>0 & p2==0){
          gen<-append(gen,'0')
        }else if (p2>0 & p1==0){gen<-append(gen,'1')}
        else{gen<-append(gen,NA)}
      }
    }
  }
  parents_good[i,5:ncol(parents_good)]<-gen
}

chr8 <-parents_good
mydim8 <- dim(chr8)
#save.image('parents_good_Asue_scaffold2.RData')
########################################################################################################################

########After doing each chromosome run HMM#######################

library(HMM)


### load SNP data

load("parents_good_Asue_scaffold13.RData")
mychr_snps <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13) 
colnames(mychr_snps)[1] <- "CHROM"
colnames(mychr_snps)[2] <- "POS"
mychr_snps <- mychr_snps[,c(1,2,5:ncol(mychr_snps))]
#Remove columns more than 90% NAs = 185 individuals left
mychr_snps <- mychr_snps[, -which(colMeans(is.na(mychr_snps)) > 0.90)]

save.image('parents_good_Asue_scaffold13.RData')


mysnpnames <- paste(mychr_snps$CHROM, mychr_snps$POS, sep="_")
rownames(mychr_snps) <- mysnpnames


my_cross2 <- apply(mychr_snps[,3:ncol(mychr_snps)],1,function(x){
  tf <- factor(x,levels=c(0,0.5,1))
  tabl1 <- table(tf)
  df <- as.data.frame(tabl1)
  return(df)
})
  
ugh <- lapply(my_cross2, function(x)  {x[,2]})
blah <- do.call(rbind, ugh)
colnames(blah) <- c("0", "0.5", "1")
blah2 <- data.frame(blah)
rownames(blah2) <- mysnpnames
blah2$ratio <- blah2$X0/(blah2$X0+blah2$X1)
blah2$called <- blah2$X0+blah2$X1
blah.keep <- blah2[blah2$ratio <=0.55 & blah2$ratio >=0.45 & blah2$called >= 20,]

mychr_snps2 <- mychr_snps[rownames(mychr_snps)%in%rownames(blah.keep),]
snps.no.te.s <- split(mychr_snps2, mychr_snps2$CHROM)


par(mfcol=c(1,13))


for (up in 1){
  up.chr <- up
  up.dat <- snps.no.te.s[[up.chr]]
  up.line <- 101
  up.snps <- up.dat[,c(1:2,up.line+2)]
  up.snps <- up.snps[is.na(up.snps[,3])==FALSE,]
  
  gt.hmm <- initHMM(c("ref", "het", "alt"), c("0","0.5","1"), transProbs=matrix(c(.9989,.001,.0001,.0001,.9998,.0001,.0001,.001,.9989),3), emissionProbs=matrix(c(.599,.4,.001,.1,.8,.1,.001,.4,.599),3))
  
  obs <- as.character(up.snps[,3])
  obs<-obs[is.na(obs)==FALSE]
  gt.vit <- viterbi(gt.hmm, obs)
  
  
  out.dat <- cbind(up.snps, gt.vit)
  par(mfcol=c(1,1))
  #plot(x=1:nrow(out.dat), out.dat[,3], type='l')
  
  gt.vit <- as.factor(gt.vit)
  #gt.vit <- as.numeric(gt.vit)
  gt.vit <- gsub("ref",0,gt.vit);gt.vit <- gsub("het",0.5,gt.vit);gt.vit <- gsub("alt",1,gt.vit)
  gt.vit <- as.numeric(gt.vit)
  plot(up.snps$POS, gt.vit, col="blue", ylim=c(-0.5,1))
  points(up.snps$POS, abs(as.numeric(up.snps[,3]))-0.1, col="red")
  legend(x="bottom", col=c("blue", "red"), legend=c("genotype","snp"), pch=19,)

hmm.gt <- function(line.no, snps){
  out.gts.s <- as.list(1)
  up.line <- line.no
  for (up.chr in 1){
    #print(up)
    up.dat <- snps.no.te.s[[up.chr]]
    up.snps <- up.dat[,c(1:2,up.line+2)]
    up.snps <- up.snps[is.na(up.snps[,3])==FALSE,]
    
    ########################## initialize HMM################################
    
    gt.hmm <- initHMM(c("ref", "het", "alt"), c("0","0.5","1"), transProbs=matrix(c(.9989,.001,.0001,.0001,.9998,.0001,.0001,.001,.9989),3), emissionProbs=matrix(c(.599,.4,.001,.1,.8,.1,.001,.4,.599),3))
    obs <- as.character(up.snps[,3])
    gt.vit <- viterbi(gt.hmm, obs)
    
    out.dat <- cbind(up.snps, gt.vit)
    out.dat2 <- cbind(data.frame(out.dat$CHROM), data.frame(out.dat$POS), gt.vit)
    out.dat2$gt.vit <- gsub("ref", 0, out.dat2$gt.vit); out.dat2$gt.vit <- gsub("het", 0.5, out.dat2$gt.vit); out.dat2$gt.vit <- gsub("alt", 1, out.dat2$gt.vit)
    names(out.dat2)[1:3] <- names(out.dat)[1:3]
    rownames(out.dat2) <- paste(out.dat2$CHROM, out.dat2$POS, sep="_")
    out.gts.s[[up.chr]] <- out.dat2
    
  }
  return(out.gts.s)
}



all.gts <- as.list(3:(ncol(snps.no.te.s[[1]])-2))
for (up.l in c(3:(ncol(snps.no.te.s[[1]])-2))){
  print(up.l)
  l.gts <- hmm.gt(line.no=up.l, snps=snps.no.te.s)
  all.gts[[up.l]] <- l.gts
}
save.image('parents_good_Asue_scaffold13.RData')

###Windows######
########### these need to be genotyped in windows
########### easy window-based filtering
########### variables : window size, prop of SNPs that need to be one gt to call window, min #SNPs in window
#Genomic windows Asue_genome.180819.HiCcor.masked.fasta.fai
window.size <- 500000 # also try 1Mb
#min.snp.no <- 10 #try 20 SNPs
#window.size = 1000000
min.snp.no = 10
#prop.gt <- .50
#window.size = 250000
#min.snp.no = 6
prop.gt <- .50

index <- read.table("/Users/robin.burns/Documents/002Asuecica/004WGA/WGA/Asue_genome.HiCGeneticMap.270919.fasta.fai")
wins <- as.list(1)
for(up in 1){
  up.m <- 19280462#index[up,2]
  up.start <- seq(1,up.m,window.size)
  up.stop <- up.start+(window.size-1)
  up.stop[length(up.stop)]<-up.m
  up.win <- cbind(up.start,up.stop)
  wins[[up]] <- up.win
}


#### step two: go through each sample, each chromosome
#### each window, call gt of window

line.gts <- c(3:length(all.gts))

line.gts.out <- lapply(wins, function(x){
  dat.m <- matrix(data=NA,nrow=nrow(x), ncol=ncol(mychr_snps2))
  return(dat.m)
})

for(up.s in line.gts){
  print(up.s)
  up.gts <- all.gts[[up.s]]
  for(up.c in 1){  
    print(up.c)
    uu.g <- up.gts[[up.c]]
    uu.w <- wins[[up.c]]
    all.w.gts <- rep(NA,nrow(uu.w))
    for(up.w in 1:nrow(uu.w)){
      up.pos <- seq(uu.w[up.w,1], uu.w[up.w,2])
      uuu.g <- uu.g[uu.g$POS%in%up.pos,]
      if(nrow(uuu.g)<=min.snp.no){w.gt <- NA}
      else{
        uuu.table <- table(uuu.g[,3])
        uuu.table <- uuu.table[order(uuu.table,decreasing=TRUE)]
        win.total <- sum(uuu.table)
        win.prop <- uuu.table[1]/win.total
        if(win.prop>=prop.gt){w.gt <- names(win.prop)} else {w.gt <- NA}
      }
      all.w.gts[up.w]<-w.gt
    }
    line.gts.out[[up.c]][,up.s]<-all.w.gts	
  }
}


#Add names
s.names <- colnames(mychr_snps2[,1:185])
  
  test <- lapply(line.gts.out, function(x){
    colnames(x)<-s.names
    return(x)
  })

  test <- do.call(rbind, test)
  
  ### get coordinates of these windows
  win.co <- wins
  win.pos <- lapply(as.list(1),function(x){
    up.dat <- win.co[[x]]
    up.chr <- rep(x,nrow(up.dat))
    up.dat <- cbind(up.chr, up.dat)
    return(up.dat)
  })
  win.pos <- do.call(rbind,win.pos)
  colnames(win.pos)<-c("chr", "window.start","window.stop")
  
  window.data <- cbind(win.pos, test)
  window.data <- window.data[,c(1,2,3,6:ncol(window.data))]
  write.csv(window.data, file="ASS3AS150F2.HMM500kb19SNPs50gt.window2020.csv",quote=FALSE,row.names=FALSE)
  window.data <- read.csv(file='ASS3AS150F2.HMM500kb19SNPs50gt.window2020.csv', header=T, stringsAsFactors = F)
  
  windowchr1 <- window.data[window.data$chr==1,]

  my_tablechr1 <- apply(windowchr1[,4:ncol(windowchr1)],1,function(x){
    tf <- factor(x,levels=c(0,0.5,1))
    tabl1 <- table(tf)
    df <- as.data.frame(tabl1)
    return(df)
  })
  
  ugh_chr1 <- lapply(my_tablechr1, function(x)  {x[,2]})
  blah_chr1 <- do.call(rbind, ugh_chr1)
  colSums((blah_chr1))

  
  dim(windowchr1)

pdf(file='Chr8_recom_2020.pdf',height = 15,width=15)
par(mfrow=c(10,4))
for (ind in 10:50){
  mytest=data.frame(cbind(windowchr1$window.start, windowchr1[,ind]))
  #mytest$snptoplot=gsub("H", 0.5, mytest[,2])
  #mytest$snptoplot2=gsub("A", 1, mytest[,3])
  #mytest$snptoplot3=gsub("B", 0, mytest[,4])
  plot(mytest[,2]~mytest[,1], ylim=c(0,1), col='blue')
}
dev.off()


pdf(file='Chr8_freqABH_2020.pdf')
freqA<-c();freqB<-c();freqH<-c();
for (i in 1:(dim(windowchr1)[1]-3)){
  freqA<-append(freqA, length(which(windowchr1[i,]==0))/length(windowchr1[i,-is.na(windowchr1[i,])]))
  freqB<-append(freqB, length(which(windowchr1[i,]==1))/length(windowchr1[i,-is.na(windowchr1[i,])]))
  freqH<-append(freqH, length(which(windowchr1[i,]==0.5))/length(windowchr1[i,-is.na(windowchr1[i,])]))
}
plot(x=1:length(freqH),y=freqH, type='l', col='black',ylim=c(0,1))
lines(x=1:length(freqA),y=freqA, col='red')
lines(x=1:length(freqB),y=freqB, col='blue')
dev.off()


#format chromosome to work with  R/qtl
windowchr1 <- window.data[window.data$chr==1,]
my_genome=data.frame(windowchr1[,4:ncol(windowchr1)])
my_genome <- apply(my_genome, 2, function(x) gsub(0.5,"H",x, perl=T));my_genome <- apply(my_genome, 2, function(x) gsub("1","B",x, perl=T));
my_genome <- apply(my_genome, 2, function(x) gsub("0","A",x, perl=T))
my_genome <- data.frame(my_genome)
my_genome$autosome=rep("autosome", nrow(my_genome))
#my_genome$ASS3a=rep("A", nrow(my_genome))
#my_genome$AS150=rep("B", nrow(my_genome))
my_genome=t(my_genome[,c(ncol(my_genome),1:ncol(my_genome)-1)])
colnames(my_genome)<-paste(windowchr1$chr, windowchr1$window.start, sep="_")
library(data.table)
fwrite(as.data.frame(my_genome), file='chr8_good.csv', na = "NA",sep=',', col.names=T, row.names=T, quote=F)

save.image('parents_good_Asue_scaffold13.RData')

library(qtl)
mapthis <- read.cross("csv", getwd(),"chr8_good.csv", na.strings="NA",estimate.map=FALSE)

summary(mapthis)



##Drop markers with a lot of missing data
#To omit the markers with lots of missing data, we first need to identify the names of the
#markers. We then use drop.markers().
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 80])
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)
pdf("Chr8_XOcounts_2020.pdf", height=7, width=12)
par(mfrow=c(1,2), las=1)
h6=hist(countXO(mapthis), breaks=183, main="XO events for chr8")
plotMissing(mapthis)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
mapthis <- fill.geno(mapthis) #increases crossover events do not do imputation
hist(countXO(mapthis), breaks=183, main="XO events for chr8 after imputation")
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
plotMissing(mapthis)
summary(mapthis)
#par(mfrow=c(1,2), las=1)
dev.off()



mapthis <- est.rf(mapthis)


rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
gt <- geno.table(mapthis)

pdf("Chr8_LODRF_2020.pdf", height=7, width=8)

hist(log(gt$P.value))

gt<-geno.table(mapthis)
plot(x=c(1:dim(gt)[1]), log(gt$P.value))
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=3)
table(lg[,2])
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=3, reorgMarkers=TRUE)
plotRF(mapthis, alternate.chrid=TRUE)
dev.off()



mapthis1 <- orderMarkers(mapthis, chr=1) #need to try multiple as there is stochasticisty to the marker order

plotRF(mapthis1, alternate.chrid=TRUE)

pull.map(mapthis1, chr=1)
plotMap(mapthis1, show.marker.names=T)



compareorder(mapthis1, chr=1, c(1:3,8,7,11,14,9,13,10,12,4,5,6,15,16:25), error.prob=0.02) #a lot of manual work to pick the best map 
#look carefully at plotRF and manually change orders as these are the diagnostic plots, straight yellow line
mapthis2=switch.order(mapthis1, chr=1,c(1:3,8,7,11,14,9,13,10,12,4,5,6,15,16:25), error.prob=0.02) #flip it around to match assembly

plotMap(mapthis2, show.marker.names=T)
plotRF(mapthis2, alternate.chrid=TRUE)
mychr2_markers=pull.map(mapthis2, chr=1)
mapthis2=mapthis1 #If it is correct at first

pdf("Chr8_GeneticMap_OrderedmarkersLOFRF_19indcalls_imputeafterfilt_2020.pdf", height=7, width=8)
plotRF(mapthis2, alternate.chrid=TRUE)
pull.map(mapthis2, chr=1)
plotMap(mapthis2, show.marker.names=T)
dev.off()



save.image('parents_good_Asue_scaffold13.RData')

library(qtl)
gt <- geno.table(mapthis)
gt_new <- geno.table(mapthis2)

library(genoPlotR)
#start1<-unlist(strsplit(rownames(gt), split='_'))
start1 <- gsub("1_","", rownames(gt))
#start1<-as.numeric(start1[seq_along(start1) %% 3 ==2])
start1<-as.numeric(start1)
end1<-start1+20000
start2<-gsub("9_","", rownames(gt_new))
start2<-as.numeric(start2)
end2<-start2+20000
comparison2<-as.data.frame(cbind(start1,end1,start2,end2))
comparison2$direction<-rep(1,dim(comparison2)[1])
comparison2$col<-rep(c('black'),dim(comparison2)[1])
#comparison2 <- comparison2[c(37:1),]
comparison2<-as.comparison(comparison2)
####plot rearrengements
chr="Asue_scaffold9"
data<-read.csv(paste('~/Desktop/Asue_AtAa_170919_dnadiff1kb.coords.9'), header=F, sep=',')
data<-data[data$V5=="scaffold_4",]
colnames(data)<-c('ref_start','ref_end','query_start','query_end','ref')
#data$query_start=19000000-data$query_start
#data$query_end=19000000-data$query_end
data$ref_start <- as.numeric(as.character(data$ref_start))
names <- c("chr")
starts <- c(1)
ends <- c(max(data$query_end))
strands <- c("+")
cols <- c("grey")
dna_seg1 <- dna_seg(data.frame(name=names, start=starts, end=ends,strand=strands, col=cols))
dna_seg2 <- dna_seg(data.frame(name=names, start=starts, end=ends,strand=strands, col=cols))
dna_seg3 <- dna_seg(data.frame(name=names, start=starts, end=ends,strand=strands, col=cols))
comparison1 <- as.comparison(data.frame(start1=data$ref_start, end1=data$ref_end,start2=data$query_start, end2=data$query_end))

comparison1$col<-'grey'
comparison1[which(comparison1$end2<comparison1$start2),'col']<-'#B22222'
pdf("Chr9_genoplotR_finalmap.pdf", height=10, width=14)
plot_gene_map(dna_segs=list(dna_seg1, dna_seg2,dna_seg3),comparisons=list(comparison1, comparison2),scale = F)
sum_as <- comparison1$end1-comparison1$start1
hist(abs(sum_as), main='hist of alignment lengths')
#save.image('parents_good_Asue_scaffold13.RData')
dev.off()


cmmb=read.table(file='ChrAll_CMvsMB.txt', header=F)



cmmb_a=cmmb[cmmb$V1 %in% c("chr1","chr4","chr9","chr11", "chr12", "chr13"),]

cmmb_b=cmmb[!(cmmb$V1 %in% c("chr1","chr4","chr9","chr11", "chr12", "chr13")),]

cmmb_a_1=cmmb_a[cmmb_a$V1 %in% c("chr1"),]
cmmb_a_1$V2=rev(cmmb_a_1$V2)
cmmb_a_4=cmmb_a[cmmb_a$V1 %in% c("chr4"),]
cmmb_a_4$V2=rev(cmmb_a_4$V2)
cmmb_a_9=cmmb_a[cmmb_a$V1 %in% c("chr9"),]
cmmb_a_9$V2=rev(cmmb_a_9$V2)
cmmb_a_11=cmmb_a[cmmb_a$V1 %in% c("chr11"),]
cmmb_a_11$V2=rev(cmmb_a_11$V2)
cmmb_a_12=cmmb_a[cmmb_a$V1 %in% c("chr12"),]
cmmb_a_12$V2=rev(cmmb_a_12$V2)
cmmb_a_13=cmmb_a[cmmb_a$V1 %in% c("chr13"),]
cmmb_a_13$V2=rev(cmmb_a_13$V2)

cmmb_a=rbind(cmmb_a_1,cmmb_a_4,cmmb_a_9,cmmb_a_11,cmmb_a_12,cmmb_a_13)

cmmb=rbind(cmmb_a, cmmb_b)

cmmb_t=cmmb[cmmb$V1%in%c("chr1", "chr3", "chr4", "chr5"),]

cmmb_a=cmmb[cmmb$V1%in%c("chr6", "chr7", "chr8", "chr9","chr10", "chr11", "chr12", "chr13"),]

cmmb_6=cmmb[cmmb$V1%in%c("chr6"),]



cmmb_t$V1=factor(cmmb_t$V1, levels=c("chr1", "chr3", "chr4", "chr5"))
                                   
cmmb_a$V1=factor(cmmb_a$V1, levels=c("chr6", "chr7", "chr8", "chr9","chr10","chr11","chr12", "chr13"))


gm1=ggplot(cmmb_t, aes(x=V2,y=V3, color=V1)) + geom_line(size=1.2) +
  scale_color_manual(name='Chr',values=c("#1b9e77","#d95f02","#7570b3","#e7298a")) +
  ylab('genetic distance (cM)') +
  xlab('physical distance (Mb)') +
  ggtitle(expression(italic('A. thaliana')*' sub-genome'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y = element_text(size=9),axis.text.x = element_text(size=9))


gm2=ggplot(cmmb_a, aes(x=V2,y=V3, color=V1)) + geom_line(size=1.2) +
  scale_color_manual(name='Chr',values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#fb9a99",'#e31a1c','#fdbf6f',
                              '#ff7f00')) +
  theme_bw() +
  ggtitle(expression(italic('A. arenosa')*' sub-genome'))+
  ylab('genetic distance (cM)') +
  xlab('physical distance (Mb)') +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y = element_text(size=9),axis.text.x = element_text(size=9))

library(ggpubr)
pdf("~/Desktop/GeneticMap_CMvsMB.pdf", width=7.8, height=4)
ggarrange(gm1,gm2, labels=c('a','b'), nrow=1, ncol=2, align='hv')
dev.off()

