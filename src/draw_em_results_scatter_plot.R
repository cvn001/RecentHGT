library(ggplot2)

mydata <- read.table('trimmed_strain_pair_hgts.txt', sep='\t', header=T)
p <- ggplot(mydata, aes(x=mydata$ani, y=mydata$hgt.num, color=mydata$type)) +
  geom_point(alpha=0.5, size=4) + scale_size_area() + scale_color_brewer(palette='Set1') + 
  ylab('Number of HGT') + xlab('ANIm of Strain Pair') + theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(color='Strains') + theme(axis.text.y=element_text(face='bold', size=10)) + 
  theme(axis.title.y=element_text(face='bold', size=12))+ theme(axis.text.x=element_text(face='bold', size=10)) + 
  theme(axis.title.x=element_text(face='bold', size=12)) + 
  ggtitle('Numbers of Recent HGT of All Strain Pairs') + theme(plot.title=element_text(face='bold', size=15)) +
  geom_hline(yintercept=40, linetype='dashed', size=1) + 
  theme(legend.title=element_text(face='bold', size=10), legend.text=element_text(face='bold', size=10))
ggsave('all_em_results_scatter_plot.png')