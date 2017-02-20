library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
data <- read.csv(input_file, header=TRUE, sep = '\t')
attach(data)
p <- ggplot(data, aes(x=Identity)) + 
  geom_histogram(bins=60) + 
  facet_wrap(~Pair. scales="free_y") + 
  ylab("Gene Counts") + 
  xlab("Sequence Identity")
ggsave(output_file)
