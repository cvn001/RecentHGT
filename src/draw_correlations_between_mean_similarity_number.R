library('ggplot2')
setwd("D:/Tong_Project/pan_genome/Strains_for_pairwise_alignment/strain_pair_hgt_identity_number_pSym")
mydata <- read.table('combined_results.txt', sep='\t', header=T)
attach(mydata,warn.conflicts = F)
p <- ggplot(mydata, aes(x=mean, y=num)) + 
  geom_point() + 
  stat_smooth(method=lm, color='black', se=T, level=0.95) +
  ylab("Number of HGTs") + ggtitle('Correlations between Mean Similarity and Numbers of HGTs') + 
  xlab("Mean Similarity of HGTs") +
  facet_wrap(~strain) + 
  theme(strip.text=element_text(size=rel(1.5)), strip.background=element_rect(size=rel(1))) + 
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14)) + 
  theme(plot.title=element_text(size=rel(1.6), face='bold'))
ggsave('correlations_between_mean_similarity_number.png')

   