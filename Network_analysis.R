library(igraph) # Load the igraph package

C=read.table("F:/parasite RNAseq Data Analysis/Result/Cytoscape_Net.txt")

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
# V(net)$community <- cfg$membership
# colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
# plot(net, vertex.color=colrs[V(net)$community])

ceb <- cluster_edge_betweenness(net)
dendPlot(ceb, mode="hclust")
plot(ceb, net)

## neighbor

neighbors(net, "Ccl4", mode="all")


###  Plot the temporal profile of larger degree nodes
# Tyrobp Ctss B2m Fermt3
dev.off()
plot(seq(0,21,1),Node_gene["Tyrobp",2:23]*max(DEG["Tyrobp",]),type='l',lty=2,main="Tyrobp",sub='',xlim=c(0,21), ylim=c(min(DEG["Tyrobp",]), max(DEG["Tyrobp",])),xlab="Time (Days)",ylab='Expression',xaxt="n",cex.axis=1.5,cex.lab=2,cex.main=2.5)
par(new=TRUE)
plot(c(0,2,7,14,21),DEG["Tyrobp",],type='p',lwd=2,main="Tyrobp",sub='',xlim=c(0,21), ylim=c(min(DEG["Tyrobp",]), max(DEG["Tyrobp",])),xlab="Time (Days)",ylab='Expression',xaxt="n",cex.axis=1.5,cex.lab=2,cex.main=2.5)
axis(1,at=c(0,2,7,14,21),las=1,cex.axis=1.5,cex.lab=2)

dev.off()
plot(seq(0,21,1),Node_gene["Ctss",2:23]*max(DEG["Ctss",]),type='l',lty=2,main="Ctss",sub='',xlim=c(0,21), ylim=c(min(DEG["Ctss",]), max(DEG["Ctss",])),xlab="Time (Days)",ylab='Expression',xaxt="n",cex.axis=1.5,cex.lab=2,cex.main=2.5)
par(new=TRUE)
plot(c(0,2,7,14,21),DEG["Ctss",],type='p',lwd=2,main="Ctss",sub='',xlim=c(0,21), ylim=c(min(DEG["Ctss",]), max(DEG["Ctss",])),xlab="Time (Days)",ylab='Expression',xaxt="n",cex.axis=1.5,cex.lab=2,cex.main=2.5)
axis(1,at=c(0,2,7,14,21),las=1,cex.axis=1.5,cex.lab=2)


dev.off()
plot(seq(0,21,1),Node_gene["B2m",2:23]*max(DEG["B2m",]),type='l',lty=2,main="B2m",sub='',xlim=c(0,21), ylim=c(min(DEG["B2m",]), max(DEG["B2m",])),xlab="Time (Days)",ylab='Expression',xaxt="n",cex.axis=1.5,cex.lab=2,cex.main=2.5)
par(new=TRUE)
plot(c(0,2,7,14,21),DEG["B2m",],type='p',lwd=2,main="B2m",sub='',xlim=c(0,21), ylim=c(min(DEG["B2m",]), max(DEG["B2m",])),xlab="Time (Days)",ylab='Expression',xaxt="n",cex.axis=1.5,cex.lab=2,cex.main=2.5)
axis(1,at=c(0,2,7,14,21),las=1,cex.axis=1.5,cex.lab=2)


dev.off()
plot(seq(0,21,1),Node_gene["Fermt3",2:23]*max(DEG["Fermt3",]),type='l',lty=2,main="Fermt3",sub='',xlim=c(0,21), ylim=c(min(DEG["Fermt3",]), max(DEG["Fermt3",])),xlab="Time (Days)",ylab='Expression',xaxt="n",cex.axis=1.5,cex.lab=2,cex.main=2.5)
par(new=TRUE)
plot(c(0,2,7,14,21),DEG["Fermt3",],type='p',lwd=2,main="Fermt3",sub='',xlim=c(0,21), ylim=c(min(DEG["Fermt3",]), max(DEG["Fermt3",])),xlab="Time (Days)",ylab='Expression',xaxt="n",cex.axis=1.5,cex.lab=2,cex.main=2.5)
axis(1,at=c(0,2,7,14,21),las=1,cex.axis=1.5,cex.lab=2)


#######  Enrichment
Gene_BP=read.table('F:/parasite RNAseq Data Analysis/Result/BP enrichment.csv',fill= T,sep=',')
Gene_BP_class=list(10)
for (i in 1:10)
{
  # Gene_BP_class[[i]]=Gene_BP[i,1]
  Gene_BP_class[[i]]=unlist(strsplit(as.character(Gene_BP[i,2]),split=',')) 
}

BP=list(10)
BP[[1]]= unique(c(unlist(Gene_BP_class[[1]]),unlist(Gene_BP_class[[2]]),unlist(Gene_BP_class[[4]]),unlist(Gene_BP_class[[5]]),unlist(Gene_BP_class[[7]])))
BP[[2]]= unique(c(unlist(Gene_BP_class[[3]])))
BP[[3]]= unique(c(unlist(Gene_BP_class[[6]]),unlist(Gene_BP_class[[8]])))
BP[[4]]= unique(c(unlist(Gene_BP_class[[9]])))
BP[[5]]= unique(c(unlist(Gene_BP_class[[10]])))


intersect(unlist(BP[[3]]),intersect(unlist(BP[[1]]),unlist(BP[[2]])))
Reduce(intersect,  list(unlist(BP[[1]]),
                        unlist(BP[[2]]),
                        unlist(BP[[3]]),
                        unlist(BP[[4]]),
                        unlist(BP[[5]])))

Union_genes=Reduce(union,  list(unlist(BP[[1]]),
                        unlist(BP[[2]]),
                        unlist(BP[[3]]),
                        unlist(BP[[4]]),
                        unlist(BP[[5]])))

Union_genes_A=Reduce(union,  list(unlist(BP[[2]]),
                                unlist(BP[[3]]),
                                unlist(BP[[4]]),
                                unlist(BP[[5]])))

setdiff(unlist(BP[[1]]),Union_genes_A)

Union_genes_D=Reduce(union,  list(unlist(BP[[1]]),
                                  unlist(BP[[2]]),
                                  unlist(BP[[3]]),
                                  unlist(BP[[5]])))

setdiff(unlist(BP[[4]]),Union_genes_D)

# sample five-set Venn Diagram  
install.packages("VennDiagram")
library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    A = unlist(BP[[1]]),
    B = unlist(BP[[2]]),
    C = unlist(BP[[3]]),
    D = unlist(BP[[4]]),
    E = unlist(BP[[5]])
  ),
  category.names = c(
    expression( bold('Immune response') ),
    expression( bold('Defence response') ),
    expression( bold('Antigen processing and presentation') ),
    expression( bold('Response to stress') ),
    expression( bold('Cell chemotaxis') )
  ),
  filename = "Venn_5set_pretty.tiff",
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = 'lzw',
  units = 'px',
  lwd = 6,
  lty = 'blank',
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  alpha = 0.50,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.5),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.2
)


# venn.diagram(list(A=1:10,B=3:8,C=6:9), fill=c("red","green","blue"), alpha=c(0.5,0.5,0.5), cex=2, filename="VennDiagram.tiff")
