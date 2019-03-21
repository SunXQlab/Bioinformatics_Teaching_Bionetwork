library(igraph) # Load the igraph package

C=read.table("Path/Cytoscape_Net.txt")

Net_R=as.data.frame(cbind(C[,1],C[,3],C[,2]))  # Input data for R analysis, data.frame


net=graph_from_data_frame(Net_R, directed = TRUE, vertices = NULL)  # 
triad_census(net) # for directed networks

#Diameter
diam <- get_diameter(net, directed=T)

# Degree
deg <- degree(net, mode="all")
# plot(net, vertex.size=deg*3)
hist(deg, breaks=1:max(deg), main="Histogram of node degree")
deg.dist <- degree_distribution(net, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

## Centrality & centralization

degree(net, mode="in")
centr_degree(net, mode="in", normalized=T)


eigen_centrality(net, directed=T, weights=NA)
centr_eigen(net, directed=T, normalized=T)

## Hubs and Authorities

hs <- hub_score(net, weights=NA)$vector
which(hs==max(hs)) # Tyrobp

as <- authority_score(net, weights=NA)$vector
which(as==max(as)) # C1qb

par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs")
plot(net, vertex.size=as*30, main="Authorities")
dev.off()


# cliques
net.sym <- as.undirected(net, mode= "each",
                         edge.attr.comb=list(weight="sum", "ignore"))

cliques(net.sym) # list of cliques
sapply(cliques(net.sym), length) # clique sizes
largest_cliques(net.sym) # cliques with max number of nodes
vcol <- rep("grey80", vcount(net.sym))
vcol[unlist(largest_cliques(net.sym))] <- "gold"
plot(as.undirected(net.sym), vertex.label=V(net.sym)$name, vertex.color=vcol)

## Community detection

cfg <- cluster_fast_greedy(as.undirected(net))
plot(cfg, as.undirected(net))
cfg$membership #label the group ID of each node


ceb <- cluster_edge_betweenness(net)
dendPlot(ceb, mode="hclust")
plot(ceb, net)

## neighbor

neighbors(net, "Ccl4", mode="all")

list(A=1:10,B=3:8,C=6:9), fill=c("red","green","blue"), alpha=c(0.5,0.5,0.5), cex=2, filename="VennDiagram.tiff")
