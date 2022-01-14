library(data.table)
library(ggplot2)
library(ggforce)
args = commandArgs(trailingOnly=TRUE)

dir<-file.path(args[1])
files<-list.files(dir, pattern=".global.dist.txt", full.names=TRUE)

all<-list()
counter<-0

for (f in files) {
  
  counter<-counter+1
  name<-basename(f)
  sample<-unlist(strsplit(name, ".", fixed=TRUE))[1]
  aligner<-unlist(strsplit(name, ".", fixed=TRUE))[2]
  tab<-fread(f, sep="\t")
  colnames(tab)<-c("chrom", "cov", "perc")
  tab$sample<-sample
  tab$aligner<-aligner
  all[[counter]]<-tab
  
}


alltab<-do.call(rbind, all)
alltab<-subset(alltab, (chrom=="total"))
ifelse(alltab$cov<=100, TRUE,FALSE)->zoom
alltab$zoom<-zoom



p<-ggplot(alltab, aes(x=cov, y=perc,col=sample,linetype=aligner)) + geom_line() + labs(x="Depth", y=expression("Fraction of bases ">=" depth")) + theme_bw() + 
  theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal")+
  facet_zoom(x=zoom)
  
ggsave(args[2])
