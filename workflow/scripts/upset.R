library(ComplexUpset)
library(data.table)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

tab<-fread(args[1], sep="\t", header=T)
samples<-unique(unlist(lapply(strsplit(colnames(tab)[c(4:ncol(tab))],".", fixed=T), function(x) unlist(x)[1])))

p<-upset(data=tab,
         intersect=samples,
         annotations = list("SV distribution"=(ggplot(mapping=aes(fill=SVTYPE)) + geom_bar(stat='count', position='fill',) + theme_classic() + scale_y_continuous(labels=scales::percent_format()) + ylab('SV distribution') +
                                                 theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.title=element_blank())),
                            "SV length"=( ggplot(mapping=aes(y=abs(SVLEN))) + geom_boxplot(na.rm=TRUE) + scale_y_continuous(trans="log10"))),
         base_annotations=list('Intersection size'=intersection_size(text_colors=c(on_background='black', on_bar='black'), text=list(size=3,angle=45,vjust=0.1,hjust=0.1)) + ylab('Intersection size')),
         width_ratio=0.1,
         keep_empty_groups=T,
         name='Combinations',
         min_degree=1
      )
ggsave(args[2], width=20, height=15)