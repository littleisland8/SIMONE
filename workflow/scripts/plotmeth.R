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

#test <- fread("lab/meth/1S.methfreq.parsed.bed")
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



#chr1 <- subset(test,(chrom=="chr1"))
#segmentate<-cpt.meanvar(chr1$y,method="PELT",penalty = "MBIC", minseglen = 5)
#indexes<-c(0, segmentate@cpts)
#segments<-rep(segmentate@param.est$mean,diff(indexes))

#seg <- data.frame(chrom=chr1$chrom, start= chr1$start, end = chr1$end, y=segments)
#segGr <- makeGRangesFromDataFrame(seg, keep.extra.columns = T)

pdf(args[2])
kp <- plotKaryotype(genome = "hg38", plot.type = 2)
#kpDataBackground(kp, data.panel = 2, col="white")
kpPoints(kp, data=Gr, data.panel = 1, col="grey50", cex=.2, r0=0, ymin=0, ymax=1)
kpAxis(kp, ymin=0, ymax=1,data.panel = 1, side = "right", cex=.2)
kpLines(kp, data = segGr, data.panel = 1, col="darkred", ymin = 0, ymax = 1)
kpHeatmap(kp, data = Gr, data.panel = 2, r0=0, r1=0.5)
dev.off()

