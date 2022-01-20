#library(data.table)
#library(karyoploteR)
#library(changepoint)
#library(tseries)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(karyoploteR))
suppressPackageStartupMessages(library(changepoint))
suppressPackageStartupMessages(library(tseries))

args = commandArgs(trailingOnly=TRUE)
file <- fread(args[1])

colnames(file) <- c("chrom", "start", "end", "y", "sample")
Gr<-makeGRangesFromDataFrame(file[,1:4], keep.extra.columns = T)

wanted <- c(paste0("chr", seq(1,22)), "chrX", "chrY")
reslist <- list()

for (chr in wanted){
  
  sub <- subset(file, (chrom==chr))
  segmentate<-cpt.meanvar(sub$y,method="PELT",penalty = "MBIC", minseglen = 5)
  indexes<-c(0, segmentate@cpts)
  segments<-rep(segmentate@param.est$mean,diff(indexes))
  
  seg <- data.frame(chrom=sub$chrom, start= sub$start, end = sub$end, y=segments)
  reslist[[chr]] <- seg
}

segmentated <- do.call(rbind,reslist)
segGr <- makeGRangesFromDataFrame(segmentated, keep.extra.columns = T)

pdf(args[2])
kp <- plotKaryotype(genome = "hg38", plot.type = 2)
kpPoints(kp, data=Gr, data.panel = 1, col="grey50", cex=.2, r0=0, ymin=0, ymax=1)
kpAxis(kp, ymin=0, ymax=1,data.panel = 1, side = "right", cex=.2)
kpLines(kp, data = segGr, data.panel = 1, col="darkred", ymin = 0, ymax = 1)
kpHeatmap(kp, data = Gr, data.panel = 2, r0=0, r1=0.5)
dev.off()

