library(ggplot2)
mydata <- read.table('trimmed_strain_pair_hgts.txt', sep = '\t', header = TRUE)
mytitle <- expression(paste("HGTs of All ", italic("Rhizobium"), " Strain Pairs"))
p <- ggplot(mydata, aes(x=mydata$hgt.num, fill=factor(mydata$type))) + 
  geom_dotplot(method='histodot', binwidth = 6.5) + geom_rug() +
  theme(axis.title.y=element_blank()) + xlab('HGT Number') + 
  scale_y_continuous(breaks = NULL) + 
  theme(legend.position=c(1,1), legend.justification=c(1,1), axis.title.x =element_text(size=14)) + 
  labs(fill='Strains') + 
  annotate('segment', x=-4, xend=30, y=0.48, yend=0.48, size=1.2, colour='grey20', arrow=arrow(ends='both', angle=90, length = unit(.2,'cm'))) + 
  annotate('text', x=14, xend=36, y=0.53, yend=0.53, colour='grey20', label='Occasional HGT', size=5) + 
  annotate('text', x=145, y=0.85, yend=0.85, colour='grey20', label='Observably Extensive HGT', size=5) + 
  annotate('segment', x=48, xend=243, y=0.80, yend = 0.80, size=1.2, colour='grey20', arrow=arrow(ends='both', angle=90, length = unit(.2,'cm'))) + 
  annotate('Segment', x=80, xend=84, y=0.121,yend=0.05, colour='grey20', size=1.1, arrow=arrow(length = unit(0.2,'cm'))) + 
  annotate('text', x=72, y=0.145, yend=0.145, label='L43~FH23', size=4, colour='grey20') +
  theme(legend.title=element_text(face='bold', size=13), legend.text=element_text(size=12)) + 
  scale_fill_brewer() + scale_color_grey() + ggtitle(mytitle) +
  theme(plot.title=element_text(face='bold', size=18)) +
  theme(plot.title = element_text(margin=margin(b=-20, unit="pt"))) +
  geom_vline(aes(xintercept=42, color='gery30'),  linetype='dashed', size=1, show.legend=F) +
  theme(axis.text.x=element_text(size=12, face='bold'))
ggsave('strain_pair_hgt_dotplot.png')