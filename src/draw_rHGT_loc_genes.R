library(ggplot2)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
mydata <- read.table(input_file, sep='\t', header=T)
strain <- unique(mydata$query.strain)
mydata$ani <- factor(mydata$ani, levels=ordered(unique(mydata$ani)))
p <- ggplot(mydata, aes(x=mydata$ani, y=mydata$number, fill=mydata$type)) +
  geom_bar(color='grey40', stat='identity', size = 0.5, position='dodge') + labs(fill='CFN42 Replicons') +
  xlab('Query Strains') + ylab('Count') +
  theme(plot.title=element_text(face='bold'), 
        axis.title.x=element_blank(), 
        axis.title.y=element_text(face='bold'), 
        axis.text.x=element_text(face='bold', angle=40, hjust=1, vjust=1, size = 10),
        legend.title=element_text(face='bold'),
        legend.text=element_text(face='bold'),
        axis.text.y=element_text(face='bold')) +
  guides(fill=guide_legend(reverse=T)) +
  scale_fill_brewer(palette='Pastel1') +
  ggtitle('Comparing rHGTs and Specific Location Genes') +
  scale_x_discrete(labels=array(strain))
ggsave(output_file)