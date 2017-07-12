library(ggplot2)
args <- commandArgs(trailingOnly = T)
# input_file <- 'pairwise_alignment_result.txt'
input_file <- args[1]
output_file <- args[2]
distribution_title <- args[3]
data <- read.csv(input_file, header=T, sep = '\t')
attach(data, warn.conflicts = F)
p <- ggplot(data, aes(x=Similarity, fill=as.factor(Category))) +
  geom_histogram(bins=30, colour='grey20') +
  ylab("Gene Counts") + 
  xlab("Sequence Similarity") +
  ggtitle(distribution_title) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill='Category') + theme(legend.position=c(0.02,0.98), legend.justification=c(0,1)) +
  theme(plot.title=element_text(size=rel(1.35))) +
  # scale_x_continuous(breaks=seq(0, 100, 2)) +
  # theme(legend.title=element_text(face='bold'), legend.text=element_text(face='bold'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values = c('white', 'grey80'))
ggsave(output_file)
