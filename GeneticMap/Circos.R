# Package circlize installation and loading into R environment
#install.packages("circlize")
library(circlize)


##############Subgenomes of A.suecica###################
cytoband = read.table("/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/sue.genome",
                      colClasses = c("character", "numeric", "numeric"), sep = "\t")


pdf("/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/Asue_genome_circos.pdf", width=6.2025, height=8.7675)

circos.clear()
circos.par("start.degree" = 168, "gap.degree"= c(rep(5,4), 12, rep(5,7), 11),canvas.ylim=c(-1,1),canvas.xlim=c(-1,1),clock.wise=T,
           track.margin=c(0,0),track.height=0.1)


# Initilizing the ideogram/cytoban
#circos.initializeWithIdeogram(cytoband, chromosome.index=)

cytoband[[1]] = factor(cytoband[[1]], levels = c("Chr1","Chr2","Chr3","Chr4", "Chr5", "Chr13","Chr12","Chr11","Chr10","Chr9","Chr8", "Chr7", "Chr6"))

circos.genomicInitialize(cytoband, sector.names = NULL,plotType = "NULL",track.height=0.1)


circos.track(ylim = c(0, 0.05), 
             bg.col = c(rep("#276419",5), rep("#752E82",8)), 
             bg.border = NA, track.height = 0.05)


circos.axis(sector.index ="Chr1",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[1,3]),
            labels = c("", "10MB","20MB"),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr2",major.tick=TRUE,
            major.at = c(0, 10000000, cytoband[2,3]),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr3",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[3,3]),
            labels = c("", "10MB","20MB"),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr4",major.tick=TRUE,
            major.at = c(0, 10000000, cytoband[4,3]),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr5",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[5,3]),
            labels = c("", "10MB","20MB"),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr6",major.tick=TRUE,
            major.at = c(cytoband[13,3],11826118, 1826118, 0),
            labels = c("", "10MB","20MB"),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr7",major.tick=TRUE,
            major.at = c(cytoband[7,3],4568019,0),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr8",major.tick=TRUE,
            major.at = c(cytoband[8,3],9280462,0),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr9",major.tick=TRUE,
            major.at = c(cytoband[9,3],7283699,0),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr10",major.tick=TRUE,
            major.at = c(cytoband[10,3],6176792,0),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr11",major.tick=TRUE,
            major.at = c(cytoband[11,3],9002944,0),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr12",major.tick=TRUE,
            major.at = c(cytoband[12,3],9590761,0),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Chr13",major.tick=TRUE,
            major.at = c(cytoband[13,3],5806088,0),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))


#Density

#mygenes=read.table(file='/Users/robin.burns/Documents/002Asuecica/004WGA/Circos/augustus_hints3sourcesAs_AtAa_newBusco.filtAll.gene.bed', header=F)
#circos.genomicDensity(mygenes,window.size = 1e5, col = c("#98cdb6"), track.height = 0.1)
#myTEs=read.table(file='/Users/robin.burns/Documents/002Asuecica/004WGA/Circos/Asuecica_TE_regions.TEs.bed', header=F)
#circos.genomicDensity(myTEs,window.size = 1e5, col = c("#88a1e7"), track.height = 0.1)
#mylinks=read.table(file='/Users/robin.burns/Documents/002Asuecica/004WGA/Circos/AsueAT_AsueAA_syntenyblocks_flipped345679.renamed.txt', header=T)
mygenes=read.table("/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/Asue_genes.bed")
myTEs=read.table("/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/Asue_genome_210620.full.RM_TE.bed")
mylinks=read.table("/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/Asbetween_flipped345679.txt", header=T)

rev_x = function(x, xrange = CELL_META$xlim) {
  xrange[2] - x + xrange[1]
}

inv <- "Chr6"
myTEs.rev = myTEs
l = myTEs.rev$V1 == inv; myTEs.rev$V2=as.numeric(myTEs.rev$V2); myTEs.rev$V3=as.numeric(myTEs.rev$V3);
myTEs.rev$V3[l] = rev_x(myTEs$V2[l], get.cell.meta.data("xlim",sector.index = inv))

myTEs.rev$V2[l] = rev_x(myTEs$V3[l], get.cell.meta.data("xlim",sector.index = inv))


inv <- "Chr7"
myTEs.rev = myTEs.rev
l = myTEs.rev$V1 == inv; myTEs.rev$V2=as.numeric(myTEs.rev$V2); myTEs.rev$V3=as.numeric(myTEs.rev$V3);
myTEs.rev$V3[l] = rev_x(myTEs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myTEs.rev$V2[l] = rev_x(myTEs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr8"
myTEs.rev = myTEs.rev
l = myTEs.rev$V1 == inv; myTEs.rev$V2=as.numeric(myTEs.rev$V2); myTEs.rev$V3=as.numeric(myTEs.rev$V3);
myTEs.rev$V3[l] = rev_x(myTEs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myTEs.rev$V2[l] = rev_x(myTEs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr9"
myTEs.rev = myTEs.rev
l = myTEs.rev$V1 == inv; myTEs.rev$V2=as.numeric(myTEs.rev$V2); myTEs.rev$V3=as.numeric(myTEs.rev$V3);
myTEs.rev$V3[l] = rev_x(myTEs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myTEs.rev$V2[l] = rev_x(myTEs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr10"
myTEs.rev = myTEs.rev
l = myTEs.rev$V1 == inv; myTEs.rev$V2=as.numeric(myTEs.rev$V2); myTEs.rev$V3=as.numeric(myTEs.rev$V3);
myTEs.rev$V3[l] = rev_x(myTEs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myTEs.rev$V2[l] = rev_x(myTEs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr11"
myTEs.rev = myTEs.rev
l = myTEs.rev$V1 == inv; myTEs.rev$V2=as.numeric(myTEs.rev$V2); myTEs.rev$V3=as.numeric(myTEs.rev$V3);
myTEs.rev$V3[l] = rev_x(myTEs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myTEs.rev$V2[l] = rev_x(myTEs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr12"
myTEs.rev = myTEs.rev
l = myTEs.rev$V1 == inv; myTEs.rev$V2=as.numeric(myTEs.rev$V2); myTEs.rev$V3=as.numeric(myTEs.rev$V3);
myTEs.rev$V3[l] = rev_x(myTEs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myTEs.rev$V2[l] = rev_x(myTEs$V3[l], get.cell.meta.data("xlim",sector.index = inv))


inv <- "Chr13"
myTEs.rev = myTEs.rev
l = myTEs.rev$V1 == inv; myTEs.rev$V2=as.numeric(myTEs.rev$V2); myTEs.rev$V3=as.numeric(myTEs.rev$V3);
myTEs.rev$V3[l] = rev_x(myTEs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myTEs.rev$V2[l] = rev_x(myTEs$V3[l], get.cell.meta.data("xlim",sector.index = inv))




inv <- "Chr6"
mygenes.rev = mygenes
l = mygenes.rev$V1 == inv; mygenes.rev$V2=as.numeric(mygenes.rev$V2); mygenes.rev$V3=as.numeric(mygenes.rev$V3);
mygenes.rev$V3[l] = rev_x(mygenes$V2[l], get.cell.meta.data("xlim",sector.index = inv))

mygenes.rev$V2[l] = rev_x(mygenes$V3[l], get.cell.meta.data("xlim",sector.index = inv))


inv <- "Chr7"
mygenes.rev = mygenes.rev
l = mygenes.rev$V1 == inv; mygenes.rev$V2=as.numeric(mygenes.rev$V2); mygenes.rev$V3=as.numeric(mygenes.rev$V3);
mygenes.rev$V3[l] = rev_x(mygenes$V2[l], get.cell.meta.data("xlim",sector.index = inv))
mygenes.rev$V2[l] = rev_x(mygenes$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr8"
mygenes.rev = mygenes.rev
l = mygenes.rev$V1 == inv; mygenes.rev$V2=as.numeric(mygenes.rev$V2); mygenes.rev$V3=as.numeric(mygenes.rev$V3);
mygenes.rev$V3[l] = rev_x(mygenes$V2[l], get.cell.meta.data("xlim",sector.index = inv))
mygenes.rev$V2[l] = rev_x(mygenes$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr9"
mygenes.rev = mygenes.rev
l = mygenes.rev$V1 == inv; mygenes.rev$V2=as.numeric(mygenes.rev$V2); mygenes.rev$V3=as.numeric(mygenes.rev$V3);
mygenes.rev$V3[l] = rev_x(mygenes$V2[l], get.cell.meta.data("xlim",sector.index = inv))
mygenes.rev$V2[l] = rev_x(mygenes$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr10"
mygenes.rev = mygenes.rev
l = mygenes.rev$V1 == inv; mygenes.rev$V2=as.numeric(mygenes.rev$V2); mygenes.rev$V3=as.numeric(mygenes.rev$V3);
mygenes.rev$V3[l] = rev_x(mygenes$V2[l], get.cell.meta.data("xlim",sector.index = inv))
mygenes.rev$V2[l] = rev_x(mygenes$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr11"
mygenes.rev = mygenes.rev
l = mygenes.rev$V1 == inv; mygenes.rev$V2=as.numeric(mygenes.rev$V2); mygenes.rev$V3=as.numeric(mygenes.rev$V3);
mygenes.rev$V3[l] = rev_x(mygenes$V2[l], get.cell.meta.data("xlim",sector.index = inv))
mygenes.rev$V2[l] = rev_x(mygenes$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr12"
mygenes.rev = mygenes.rev
l = mygenes.rev$V1 == inv; mygenes.rev$V2=as.numeric(mygenes.rev$V2); mygenes.rev$V3=as.numeric(mygenes.rev$V3);
mygenes.rev$V3[l] = rev_x(mygenes$V2[l], get.cell.meta.data("xlim",sector.index = inv))
mygenes.rev$V2[l] = rev_x(mygenes$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr13"
mygenes.rev = mygenes.rev
l = mygenes.rev$V1 == inv; mygenes.rev$V2=as.numeric(mygenes.rev$V2); mygenes.rev$V3=as.numeric(mygenes.rev$V3);
mygenes.rev$V3[l] = rev_x(mygenes$V2[l], get.cell.meta.data("xlim",sector.index = inv))
mygenes.rev$V2[l] = rev_x(mygenes$V3[l], get.cell.meta.data("xlim",sector.index = inv))


inv <- "Chr6"
mylinks.rev = mylinks
l = mylinks.rev$Aasue_chr == inv
mylinks.rev$Aasue_end[l] = abs(mylinks.rev$Aasue_end[l]-21826118)
mylinks.rev$Aasue_start[l] = abs(mylinks.rev$Aasue_start[l]-21826118)

inv <- "Chr7"
mylinks.rev = mylinks.rev
l = mylinks.rev$Aasue_chr == inv
mylinks.rev$Aasue_end[l] = abs(mylinks.rev$Aasue_end[l]-14568019)
mylinks.rev$Aasue_start[l] = abs(mylinks.rev$Aasue_start[l]-14568019)
  
  
inv <- "Chr8"
mylinks.rev = mylinks.rev
l = mylinks.rev$Aasue_chr == inv
mylinks.rev$Aasue_end[l] = abs(mylinks.rev$Aasue_end[l]-19280462)
mylinks.rev$Aasue_start[l] = abs(mylinks.rev$Aasue_start[l]-19280462)

inv <- "Chr9"
mylinks.rev = mylinks.rev
l = mylinks.rev$Aasue_chr == inv
mylinks.rev$Aasue_end[l] = abs(mylinks.rev$Aasue_end[l]-17283699)
mylinks.rev$Aasue_start[l] = abs(mylinks.rev$Aasue_end[l]-17283699)

inv <- "Chr10"
mylinks.rev = mylinks.rev
l = mylinks.rev$Aasue_chr == inv 
mylinks.rev$Aasue_end[l] = rev_x(mylinks$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
mylinks.rev$Aasue_start[l] = rev_x(mylinks$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr11"
mylinks.rev = mylinks.rev
l= mylinks.rev$Aasue_chr ==inv
mylinks.rev$Aasue_end[l] = rev_x(mylinks$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
mylinks.rev$Aasue_start[l] = rev_x(mylinks$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr12"
mylinks.rev = mylinks.rev
l = mylinks.rev$Aasue_chr == inv
mylinks.rev$Aasue_end[l] = rev_x(mylinks$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
mylinks.rev$Aasue_start[l] = rev_x(mylinks$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "Chr13"
mylinks.rev = mylinks.rev
l = mylinks.rev$Aasue_chr == inv
mylinks.rev$Aasue_end[l] = rev_x(mylinks$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
mylinks.rev$Aasue_start[l] = rev_x(mylinks$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))


#Chip data 
#library(data.table)
#ASS3_H3K9me2=fread(file='/Users/robin.burns/Documents/002Asuecica/006CHIP/ASS3_H3K9me2_rep1and2_Log2.bigWig.mean.bedGraph', header=F)

#inv <- "chr6"
#ASS3_H3K9me2.rev = ASS3_H3K9me2
#l = ASS3_H3K9me2.rev$Aasue_chr == inv
#ASS3_H3K9me2.rev$Aasue_end[l] = rev_x(ASS3_H3K9me2$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
#ASS3_H3K9me2.rev$Aasue_start[l] = rev_x(ASS3_H3K9me2$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

#inv <- "chr7"
#ASS3_H3K9me2.rev = ASS3_H3K9me2
#l = ASS3_H3K9me2.rev$Aasue_chr == inv
#ASS3_H3K9me2.rev$Aasue_end[l] = rev_x(ASS3_H3K9me2$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
#ASS3_H3K9me2.rev$Aasue_start[l] = rev_x(ASS3_H3K9me2$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

#inv <- "chr8"
#ASS3_H3K9me2.rev = ASS3_H3K9me2
#l = ASS3_H3K9me2.rev$Aasue_chr == inv
#ASS3_H3K9me2.rev$Aasue_end[l] = rev_x(ASS3_H3K9me2$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
#ASS3_H3K9me2.rev$Aasue_start[l] = rev_x(ASS3_H3K9me2$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

#inv <- "chr9"
#ASS3_H3K9me2.rev = ASS3_H3K9me2
#l = ASS3_H3K9me2.rev$Aasue_chr == inv
#ASS3_H3K9me2.rev$Aasue_end[l] = rev_x(ASS3_H3K9me2$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
#ASS3_H3K9me2.rev$Aasue_start[l] = rev_x(ASS3_H3K9me2$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

#inv <- "chr10"
#ASS3_H3K9me2.rev = ASS3_H3K9me2
#l = ASS3_H3K9me2.rev$Aasue_chr == inv
#ASS3_H3K9me2.rev$Aasue_end[l] = rev_x(ASS3_H3K9me2$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
#ASS3_H3K9me2.rev$Aasue_start[l] = rev_x(ASS3_H3K9me2$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

#inv <- "chr11"
#ASS3_H3K9me2.rev = ASS3_H3K9me2
#l = ASS3_H3K9me2.rev$Aasue_chr == inv
#ASS3_H3K9me2.rev$Aasue_end[l] = rev_x(ASS3_H3K9me2$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
#ASS3_H3K9me2.rev$Aasue_start[l] = rev_x(ASS3_H3K9me2$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

#inv <- "chr12"
#ASS3_H3K9me2.rev = ASS3_H3K9me2
#l = ASS3_H3K9me2.rev$Aasue_chr == inv
#ASS3_H3K9me2.rev$Aasue_end[l] = rev_x(ASS3_H3K9me2$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
#ASS3_H3K9me2.rev$Aasue_start[l] = rev_x(ASS3_H3K9me2$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))

#inv <- "chr13"
#ASS3_H3K9me2.rev = ASS3_H3K9me2
#l = ASS3_H3K9me2.rev$Aasue_chr == inv
#ASS3_H3K9me2.rev$Aasue_end[l] = rev_x(ASS3_H3K9me2$Aasue_start[l], get.cell.meta.data("xlim",sector.index = inv))
#ASS3_H3K9me2.rev$Aasue_start[l] = rev_x(ASS3_H3K9me2$Aasue_end[l], get.cell.meta.data("xlim",sector.index = inv))




circos.genomicDensity(mygenes.rev,window.size = 1e5, col = c("#98cdb6"), track.height = 0.15)
circos.genomicDensity(myTEs.rev,window.size = 1e5, col = c("#215993"), track.height = 0.15)

#circos.genomicDensity(ASS3_H3K9me2.rev,window.size = 1e5, col = c("#88a1e7"), track.height = 0.15)


#circos.link(mylinks[1,1], c(mylinks[1,2], mylinks[1,3]), mylinks[1,4], c(mylinks[1,5], mylinks[1,6]), border = 1)
mylinks1=mylinks.rev[mylinks.rev$Athsue_chr=="Chr1",]
circos.genomicLink(mylinks1[,c(2,3,4)], mylinks1[,c(1,5,6)], col = "#fdf1cb", 
                   border = NA)
mylinks2=mylinks.rev[mylinks.rev$Athsue_chr=="Chr2",]
circos.genomicLink(mylinks2[,c(2,3,4)], mylinks2[,c(1,5,6)], col = "#d6ccd8", 
                   border = NA)
mylinks3=mylinks.rev[mylinks.rev$Athsue_chr=="Chr3",]
circos.genomicLink(mylinks3[,c(2,3,4)], mylinks3[,c(1,5,6)], col = "#cdd6ad", 
                   border = NA)
mylinks4=mylinks.rev[mylinks.rev$Athsue_chr=="Chr4",]
circos.genomicLink(mylinks4[,c(2,3,4)], mylinks4[,c(1,5,6)], col = "#dbecf6", 
                   border = NA)
mylinks5=mylinks.rev[mylinks.rev$Athsue_chr=="Chr5",]
circos.genomicLink(mylinks5[,c(2,3,4)], mylinks5[,c(1,5,6)], col = "#ebb5a9", 
                   border = NA)
mylinks2=mylinks.rev[mylinks.rev$Athsue_chr=="Chr2",]
circos.genomicLink(mylinks2[,c(2,3,4)], mylinks2[,c(1,5,6)], col = "#d6ccd8", 
                   border = NA)



dev.off()
















###############Athaliana and Athaliana sub-genome of Asuecica############


cytoband = read.table("/Users/robin.burns/Documents/002Asuecica/004WGA/Circos/at_atsue.genome",
                      colClasses = c("character", "numeric", "numeric"), sep = "\t")



pdf("/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/AsueAth_circos.pdf", width=4, height=4.38375)

circos.clear()
circos.par("start.degree" = 177, "gap.degree"= c(rep(4,10)),canvas.ylim=c(-1,1),canvas.xlim=c(-1,1),clock.wise=T,
           track.margin=c(0,0),track.height=0.1)


# Initilizing the ideogram/cytoban
#circos.initializeWithIdeogram(cytoband, chromosome.index=)

cytoband[[1]] = factor(cytoband[[1]], levels = c("At1","At2","At3","At4", "At5", "As5","As4","As3","As2", "As1"))

circos.genomicInitialize(cytoband, plotType = NULL,tickLabelsStartFromZero = TRUE,
                         axis.labels.cex = 0.3*par("cex"), labels.cex = 0.9*par("cex"),track.height=0.1)
           

#circos.genomicInitialize(cytoband, sector.names = NULL,plotType = c("axis", "labels"),
#                         tickLabelsStartFromZero = TRUE,
#                         axis.labels.cex = 0.3*par("cex"), labels.cex = 0.9*par("cex"),track.height=0.1)

#"#276419"
circos.track(ylim = c(0, 0.05), 
             bg.col = c(rep("#7fbc41",5), rep("#276419",5)), 
             bg.border = NA, track.height = 0.05)


circos.axis(sector.index ="At1",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000,30000000, cytoband[6,3]),
            labels = c("", "10MB","20MB", "30MB"),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="At2",major.tick=TRUE,
            major.at = c(0, 10000000, cytoband[7,3]),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="At3",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[8,3]),
            labels = c("", "10MB","20MB"),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="At4",major.tick=TRUE,
            major.at = c(0, 10000000, cytoband[9,3]),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="At5",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[10,3]),
            labels = c("", "10MB","20MB"),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As1",major.tick=TRUE,
            major.at = c(cytoband[1,3],19820293, 9820293, 0),
            labels = c("", "10MB","20MB"),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As2",major.tick=TRUE,
            major.at = c(cytoband[2,3],9960976,0),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As3",major.tick=TRUE,
            major.at = c(cytoband[3,3],13973311,3973311,0),
            labels = c("", "10MB","20MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As4",major.tick=TRUE,
            major.at = c(cytoband[4,3],9168148,0),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As5",major.tick=TRUE,
            major.at = c(cytoband[5,3],16106147,6106147,0),
            labels = c("", "10MB","20MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))


myrelations=read.table(file='/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/Tair10_AsueAT_synteny.flipped.345.txt', header=T)
  
relations=myrelations[,c(1,5,6,2,3,4)]
colnames(relations)=c("V1", "V2", "V3", "V4", "V5", "V6")

inv <- "As1"
rev_x = function(x, xrange = CELL_META$xlim) {
  as.numeric(xrange[2]) - x + as.numeric(xrange[1])
}

relations.rev = relations
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
                                                                      
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))


inv <- "As2"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As3"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As4"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As5"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))


circos.genomicLink(relations.rev[1,1:3], relations.rev[1,4:6], col = "#e4e4e5",border = "#c0c2c4")

circos.genomicLink(relations.rev[2,1:3], relations.rev[2,4:6], col = "#e4e4e5", border = "#c0c2c4")

circos.genomicLink(relations.rev[3,1:3], relations.rev[3,4:6], col = "#e4e4e5",border = "#c0c2c4")

circos.genomicLink(relations.rev[4,1:3], relations.rev[4,4:6], col = "#e4e4e5", border = "#c0c2c4")

circos.genomicLink(relations.rev[5,1:3], relations.rev[5,4:6], col = "#e4e4e5", border = "#c0c2c4")

circos.genomicLink(relations.rev[6,1:3], relations.rev[6,4:6], col = "#e4e4e5", border = "#c0c2c4")

circos.genomicLink(relations.rev[7,1:3], relations.rev[7,4:6], col = "#e4e4e5", border = "#c0c2c4")

circos.genomicLink(relations.rev[8,1:3], relations.rev[8,4:6], col = "#e4e4e5", border = "#c0c2c4")
circos.genomicLink(relations.rev[9,1:3], relations.rev[9,4:6], col = "#e4e4e5", border = "#c0c2c4")

circos.genomicLink(relations.rev[10,1:3], relations.rev[10,4:6], col = "#e4e4e5", border = "#c0c2c4")
circos.genomicLink(relations.rev[11,1:3], relations.rev[11,4:6], col = "#e4e4e5", border = "#c0c2c4")
circos.genomicLink(relations.rev[12,1:3], relations.rev[12,4:6], col = "#e4e4e5", border = "#c0c2c4")
circos.genomicLink(relations.rev[13,1:3], relations.rev[13,4:6], col = "#e4e4e5", border = "#c0c2c4")

dev.off()



###############Alyrata and Aarenosa sub-genome of Asuecica############


cytoband = read.table("/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/al_aasue.genome",
                      colClasses = c("character", "numeric", "numeric"), sep = "\t")



pdf("/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/AsueAly_circos.pdf", width=4, height=4.38375)

circos.clear()
circos.par("start.degree" = 184, "gap.degree"= c(rep(4,7),10,rep(4,7),12),canvas.ylim=c(-1,1),canvas.xlim=c(-1,1),clock.wise=T,
           track.margin=c(0,0),track.height=0.1)


cytoband[[1]] = factor(cytoband[[1]], levels = c("Al1","Al2","Al3","Al4", "Al5", "Al6","Al7","Al8","As13", 
                                                 "As12","As11","As10","As9","As8","As7","As6"))

circos.genomicInitialize(cytoband, plotType = NULL,tickLabelsStartFromZero = TRUE,
                         axis.labels.cex = 0.3*par("cex"), labels.cex = 0.9*par("cex"),track.height=0.1)

##As=#762a83',Aa='#807dba'
circos.track(ylim = c(0, 0.05), 
             bg.col = c(rep("#807dba",8), rep("#762a83",8)), 
             bg.border = NA, track.height = 0.05)


circos.axis(sector.index ="Al1",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000,30000000, cytoband[9,3]),
            labels = c("", "10MB","20MB", "30MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Al2",major.tick=TRUE,
            major.at = c(0, 10000000, cytoband[10,3]),
            labels = c("", "10MB",""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Al3",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[11,3]),
            labels = c("", "10MB","20MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Al4",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[12,3]),
            labels = c("", "10MB","20MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Al5",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[13,3]),
            labels = c("", "10MB","20MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Al6",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[14,3]),
            labels = c("", "10MB","20MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Al7",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[15,3]),
            labels = c("", "10MB","20MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="Al8",major.tick=TRUE,
            major.at = c(0, 10000000, 20000000, cytoband[16,3]),
            labels = c("", "10MB","20MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As6",major.tick=TRUE,
            major.at = c(cytoband[1,3],11826118,1826118,0),
            labels = c("", "10MB","20MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As7",major.tick=TRUE,
            major.at = c(cytoband[2,3],4568019,0),
            labels = c("", "10MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As8",major.tick=TRUE,
            major.at = c(cytoband[3,3],4664457,0),
            labels = c("", "10MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As9",major.tick=TRUE,
            major.at = c(cytoband[4,3],7283699,0),
            labels = c("", "10MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As10",major.tick=TRUE,
            major.at = c(cytoband[5,3],6176792,0),
            labels = c("", "10MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))


circos.axis(sector.index ="As11",major.tick=TRUE,
            major.at = c(cytoband[6,3],9002944,0),
            labels = c("", "10MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))

circos.axis(sector.index ="As12",major.tick=TRUE,
            major.at = c(cytoband[7,3],9590761,0),
            labels = c("", "10MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))


circos.axis(sector.index ="As13",major.tick=TRUE,
            major.at = c(cytoband[8,3],5806088,0),
            labels = c("", "10MB", ""),minor.ticks = 0,labels.facing = "inside",labels.cex = 0.6*par("cex"))




myrelations=read.table(file='/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/Alyrata_AsueAA.syntenyblocks.flipped.124.txt', header=T)
relations=myrelations
relations=relations[,c(2,3,4,1,5,6)]
colnames(relations)=c("V4", "V5", "V6", "V1", "V2", "V3")


rev_x = function(x, xrange = CELL_META$xlim) {
  as.numeric(xrange[2]) - x + as.numeric(xrange[1])
}
inv <- "As6"
relations.rev = relations
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))

relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))


inv <- "As7"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As8"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As9"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As10"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As11"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As12"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As13"
relations.rev = relations.rev
l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))

myNs=read.table(file='/Users/robin.burns/Documents/002Asuecica/004WGA/WGAJune2020/Aly_AsueAA_N.bed', header=F)

inv <- "As6"
myNs.rev = myNs
l = myNs.rev$V1 == inv; myNs.rev$V2=as.numeric(myNs.rev$V2); myNs.rev$V3=as.numeric(myNs.rev$V3);
myNs.rev$V3[l] = rev_x(myNs$V2[l], get.cell.meta.data("xlim",sector.index = inv))

myNs.rev$V2[l] = rev_x(myNs$V3[l], get.cell.meta.data("xlim",sector.index = inv))


inv <- "As7"
myNs.rev = myNs.rev
l = myNs.rev$V1 == inv; myNs.rev$V2=as.numeric(myNs.rev$V2); myNs.rev$V3=as.numeric(myNs.rev$V3);
myNs.rev$V3[l] = rev_x(myNs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myNs.rev$V2[l] = rev_x(myNs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As8"
myNs.rev = myNs.rev
l = myNs.rev$V1 == inv; myNs.rev$V2=as.numeric(myNs.rev$V2); myNs.rev$V3=as.numeric(myNs.rev$V3);
myNs.rev$V3[l] = rev_x(myNs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myNs.rev$V2[l] = rev_x(myNs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As9"
myNs.rev = myNs.rev
l = myNs.rev$V1 == inv; myNs.rev$V2=as.numeric(myNs.rev$V2); myNs.rev$V3=as.numeric(myNs.rev$V3);
myNs.rev$V3[l] = rev_x(myNs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myNs.rev$V2[l] = rev_x(myNs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As10"
myNs.rev = myNs.rev
l = myNs.rev$V1 == inv; myNs.rev$V2=as.numeric(myNs.rev$V2); myNs.rev$V3=as.numeric(myNs.rev$V3);
myNs.rev$V3[l] = rev_x(myNs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myNs.rev$V2[l] = rev_x(myNs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As11"
myNs.rev = myNs.rev
l = myNs.rev$V1 == inv; myNs.rev$V2=as.numeric(myNs.rev$V2); myNs.rev$V3=as.numeric(myNs.rev$V3);
myNs.rev$V3[l] = rev_x(myNs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myNs.rev$V2[l] = rev_x(myNs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As12"
myNs.rev = myNs.rev
l = myNs.rev$V1 == inv; myNs.rev$V2=as.numeric(myNs.rev$V2); myNs.rev$V3=as.numeric(myNs.rev$V3);
myNs.rev$V3[l] = rev_x(myNs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myNs.rev$V2[l] = rev_x(myNs$V3[l], get.cell.meta.data("xlim",sector.index = inv))

inv <- "As13"
myNs.rev = myNs.rev
l = myNs.rev$V1 == inv; myNs.rev$V2=as.numeric(myNs.rev$V2); myNs.rev$V3=as.numeric(myNs.rev$V3);
myNs.rev$V3[l] = rev_x(myNs$V2[l], get.cell.meta.data("xlim",sector.index = inv))
myNs.rev$V2[l] = rev_x(myNs$V3[l], get.cell.meta.data("xlim",sector.index = inv))





circos.genomicDensity(myNs.rev,window.size = 1e5, col = c("#f19f32"), track.height = 0.1)

circos.genomicLink(relations.rev[1,1:3], relations.rev[1,4:6], col = "#e4e4e5",border = "#c0c2c4")
circos.genomicLink(relations.rev[2,1:3], relations.rev[2,4:6], col = "#e4e4e5",border = "#c0c2c4")
circos.genomicLink(relations.rev[4,1:3], relations.rev[4,4:6], col = "#e4e4e5",border = "#c0c2c4")
circos.genomicLink(relations.rev[6,1:3], relations.rev[6,4:6], col = "#e4e4e5",border = "#c0c2c4")
circos.genomicLink(relations.rev[8,1:3], relations.rev[8,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[9,1:3], relations.rev[9,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[10,1:3], relations.rev[10,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[11,1:3], relations.rev[11,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[11,1:3], relations.rev[11,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[12,1:3], relations.rev[12,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[14,1:3], relations.rev[14,4:6], col = "#e4e4e5",border = "#c0c2c4")  
circos.genomicLink(relations.rev[15,1:3], relations.rev[15,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[16,1:3], relations.rev[16,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[17,1:3], relations.rev[17,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[18,1:3], relations.rev[18,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[19,1:3], relations.rev[19,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[20,1:3], relations.rev[20,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[21,1:3], relations.rev[21,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[22,1:3], relations.rev[22,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[24,1:3], relations.rev[24,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[25,1:3], relations.rev[25,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[26,1:3], relations.rev[26,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[27,1:3], relations.rev[27,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[28,1:3], relations.rev[28,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[29,1:3], relations.rev[29,4:6], col = "#e4e4e5",border = "#c0c2c4")
circos.genomicLink(relations.rev[30,1:3], relations.rev[30,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[31,1:3], relations.rev[31,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[32,1:3], relations.rev[32,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[33,1:3], relations.rev[33,4:6], col = "#e4e4e5",border = "#c0c2c4")
circos.genomicLink(relations.rev[35,1:3], relations.rev[35,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[36,1:3], relations.rev[36,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[37,1:3], relations.rev[37,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[38,1:3], relations.rev[38,4:6], col = "#e4e4e5",border = "#c0c2c4") 
circos.genomicLink(relations.rev[39,1:3], relations.rev[39,4:6], col = "#e4e4e5",border = "#c0c2c4") 



circos.genomicLink(relations.rev[7,1:3], relations.rev[7,4:6], col = "#c3e0f0",border = "#64b7da") #inversion 304766
circos.genomicLink(relations.rev[23,1:3], relations.rev[23,4:6], col = "#c3e0f0",border = "#64b7da") #inversion 890781
circos.genomicLink(relations.rev[5,1:3], relations.rev[5,4:6], col = "#ebb5a9",border = "#c71f2b") #missassembly

circos.genomicLink(relations.rev[40,1:3], relations.rev[40,4:6], col = "#c3e0f0",border = "#64b7da") #inversion 806708
circos.genomicLink(relations.rev[34,1:3], relations.rev[34,4:6], col = "#c3e0f0",border = "#64b7da") #inversion 2420288
circos.genomicLink(relations.rev[3,1:3], relations.rev[3,4:6], col = "#c3e0f0",border = "#64b7da") #86953 #inversion total 4509496
circos.genomicLink(relations.rev[13,1:3], relations.rev[13,4:6], col = "#224297",border = "#224297") #translocation 
dev.off()

