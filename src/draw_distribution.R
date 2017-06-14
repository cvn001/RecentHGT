library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
distribution_title <- args[3]
data <- read.csv(input_file, header=TRUE, sep = '\t')
attach(data)
p <- ggplot(data, aes(x=Similarity)) +
  geom_histogram(bins=60) + 
  # facet_wrap(~Pair, scales="free_y") + 
  ylab("Gene Counts") + 
  xlab("Sequence Similarity") +
  ggtitle(distribution_title) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(output_file)
