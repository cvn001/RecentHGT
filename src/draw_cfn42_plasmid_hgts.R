library(ggplot2)

mydata <- read.table('CFN42_Query_Strains_HGT_Num_Statistic.txt', sep='\t', header=T)
p <- ggplot(mydata, aes(x=mydata$query, y=mydata$number, fill=mydata$plasmid)) + 
  geom_bar(color='grey40', stat='identity', size = 0.5) + labs(fill='CFN42 Plasmids') + 
  xlab('Query Strains') + ylab('Number of HGTs') + 
  theme(plot.title=element_text(face='bold'), axis.title.x=element_text(face='bold'), axis.title.y=element_text(face='bold'), axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + 
  guides(fill=guide_legend(reverse=T)) + 
  scale_fill_brewer(palette='Pastel1') + 
  ggtitle('Recent Horizontal Transfered Genes between CFN42 and Query Strains')+
  theme(axis.text.x=element_text(face='bold', size=10), axis.text.y=element_text(face='bold', size=10)) +
  theme(legend.title=element_text(face='bold'), legend.text=element_text(face='bold'))
ggsave('CFN42_Plasmids_HGTs.png')