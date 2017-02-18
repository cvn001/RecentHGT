library(igraph)
links <- read.csv("trimmed_strain_pair_hgts.txt", header=T, as.is=T)
nodes <- read.csv("strains_network_nodes.txt", header = T, as.is=T)
links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
net <- simplify(net, remove.multiple = F, remove.loops = T) 
colrs <- c("gray50", "tomato")
V(net)$color <- colrs[V(net)$strain.location]
E(net)$width <- E(net)$weight/60
par(mar=c(0,0,0,0))
plot(net, vertex.shape="none", vertex.label=V(net)$strain.name, 
       vertex.label.font=2, vertex.label.color=V(net)$color, vertex.label.cex=1, 
       edge.color="gray80", layout=layout.circle)
legend(x=1.0, y=1.0, c("Source strains","Native strains"), pch=22,
      col="#777777", pt.bg=colrs, pt.cex=1.5, cex=0.9, bty="n", ncol=1)