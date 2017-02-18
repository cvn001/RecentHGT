library(ggplot2)
library(plyr)

# rm(mydata)
mydata <- read.table('CFN42_high_identity_genes_locus_num.txt', sep='\t', header=T)
strain <- unique(mydata$query.strain)
mydata$ani <- factor(mydata$ani, levels=ordered(unique(mydata$ani)))
# ce <- ddply(mydata, "query.strain", transform, percent=mydata$num/mydata$num * 100)
p <- ggplot(mydata, aes(x=mydata$ani, y=mydata$num, fill=mydata$locus)) +
  geom_bar(color='grey40', stat='identity', size = 0.5) + labs(fill='CFN42 Replicons') +
  xlab('Query Strains') + ylab('Count') +
  theme(plot.title=element_text(face='bold'), axis.title.x=element_blank(), axis.title.y=element_text(face='bold'), axis.text.x=element_text(face='bold', angle=40, hjust=1, vjust=1, size = 10)) +
  guides(fill=guide_legend(reverse=T)) +
  scale_fill_brewer(palette='Pastel1') +
  ggtitle('The Number of Highly Similar Genes on either Chromosome or Plasmid') +
  scale_x_discrete(labels=array(strain)) + 
  theme(legend.title=element_text(face='bold'), legend.text=element_text(face='bold')) +
  theme(axis.text.y=element_text(face='bold'))
ggsave('CFN42_High_Identity_Genes_Locus.png')