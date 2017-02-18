library(ggplot2)
# args <- commandArgs(trailingOnly = TRUE)
# input_file <- args[1]
# output_file <- args[2]
mydata <- read.csv('CFN42_FH14_FH23.txt', header=TRUE, sep = '\t')
attach(mydata,warn.conflicts = F)
p <- ggplot(mydata, aes(x=Identity)) + 
     geom_histogram(bins=40) + 
     ylab("Number of Genes") + 
     facet_wrap(~Pair) + 
     xlab("Sequence Identity") +
     theme(strip.text=element_text(size=rel(1.5))) + 
     theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))
ggsave('CFN42_FH14_FH23_pSym.png')