#coverage plot
library(reshape2)
library(data.table)
library(zoo)
library(ggplot2)

plotcoverage <- function(file) {
depth <- fread(file, sep="\t", header=F)
colnames(depth)=c("Chr", "locus", "coverage")             
depth.average<-setDT(depth)[, .(
				  window.start = rollapply(locus, width=500000, by=500000, FUN=min, align="left", partial=TRUE),
				    window.end = rollapply(locus, width=500000, by=500000, FUN=max, align="left", partial=TRUE),
				      coverage = rollapply(coverage, width=500000, by=500000, FUN=mean, align="left", partial=TRUE)
				  ), .(Chr)]

depth.average$Chr=factor(depth.average$Chr, levels=c("Asue_scaffold1", "Asue_scaffold2", "Asue_scaffold3", "Asue_scaffold4", "Asue_scaffold5", "Asue_scaffold6", "Asue_scaffold7", "Asue_scaffold8", "Asue_scaffold9", "Asue_scaffold10", "Asue_scaffold11", "Asue_scaffold12", "Asue_scaffold13"))
depth.average.plot <- ggplot(depth.average, aes(x=window.end, y=coverage, colour=Chr)) + 
  #geom_point(shape = 20, size = 1) +
  geom_line(size=0.8) +
  scale_x_continuous(name="Genomic Position (bp)", labels = scales::scientific, breaks = c(0, 10000000, 20000000, 30000000)) +
  scale_y_continuous(name="Average Coverage Depth", limits=c(0, 180), breaks = c(0,30,50,70,90,120,150,180)) +
  scale_color_manual(values=c(rep("#276419",5), rep("#762a83",8))) +
    theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(legend.position="none")
  X.p5 <<- depth.average.plot + facet_grid(. ~ Chr, space="free_x", scales = "free_x")+
    theme(panel.spacing.x = grid::unit(0, "cm")) 
}


#  mypdf=paste(file, ".pdf", sep='')
#  pdf(mypdf, width=13, height=6)
#  X.p5
#  dev.off()

setwd("/scratch-cbe/users/robin.burns/021MapAs_2020/coverage_plots/genome")
args <- commandArgs(TRUE)
plotcoverage(file=args[1]) 
mypdf=args[2]
pdf(mypdf, width=10.8, height=4)
X.p5
dev.off()

