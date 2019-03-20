

##### read DEGs 
setwd("F:/parasite RNAseq Data Analysis/Code")  

install.packages("ggm")
library(ggm)

source("pcor.R")

if(!require("ppcor")){
  install.packages("ppcor", repos='http://cran.us.r-project.org')
  library(ppcor)
}

DEG=read.table("F:/parasite RNAseq Data Analysis/Result/DEGs.xls",fill= FALSE)


##### Initial PPI network
PPI_pair=read.table("F:/parasite RNAseq Data Analysis/Result/PPI_interactions.txt",header = FALSE)
NetGene12=PPI_pair[,1:2]
NetGene=union(NetGene12[,1],NetGene12[,2])
Netsize=length(NetGene) 

Node_gene=DEG[NetGene,]
  
PPI=matrix(0,nrow=Netsize,ncol=Netsize)
for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    for (k in 1:dim(NetGene12)[1])
    {
      #if (is.na(prod(which(as.matrix(PPI_pair)[k,]==c(names_SR[i],names_SR[j]))))==0)
      if (sum(c(NetGene[i],NetGene[j])==as.matrix(NetGene12)[k,])==2)
        PPI[i,j]=PPI[i,j]+1
      else
        PPI[i,j]=PPI[i,j]+0  

    }
    
  }  
}

for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    PPI[i,j]=max(PPI[i,j],PPI[j,i])
  }
}

rownames(PPI)=NetGene
colnames(PPI)=NetGene
sum(PPI)  #=415  edges
PPI["Serping1","Fn1"]
PPI["Fn1","Serping1"]
PPI["Serping1","Oasl2"]
PPI["Mmp3","Oasl2"]

setwd("F:/parasite RNAseq Data Analysis/Result")  
write.table(PPI, file="PPI_code.xls",col.names = NA,sep = "\t") 




##  Partial correlation for Sensitive gene

DEG_cluster_1_con=read.csv('F:/parasite RNAseq Data Analysis/Result/DEG_cluster_1_con.csv',fill= T,header = F)
DEG_cluster_2_con=read.csv('F:/parasite RNAseq Data Analysis/Result/DEG_cluster_2_con.csv',fill= T,header = F)
DEG_cluster_3_con=read.csv('F:/parasite RNAseq Data Analysis/Result/DEG_cluster_3_con.csv',fill= T,header = F)
DEG_cluster_4_con=read.csv('F:/parasite RNAseq Data Analysis/Result/DEG_cluster_4_con.csv',fill= T,header = F)
DEG_cluster_5_con=read.csv('F:/parasite RNAseq Data Analysis/Result/DEG_cluster_5_con.csv',fill= T,header = F)
DEG_cluster_6_con=read.csv('F:/parasite RNAseq Data Analysis/Result/DEG_cluster_6_con.csv',fill= T,header = F)
DEG_cluster_7_con=read.csv('F:/parasite RNAseq Data Analysis/Result/DEG_cluster_7_con.csv',fill= T,header = F)


DEG_cluster_con=rbind(DEG_cluster_1_con,DEG_cluster_2_con,DEG_cluster_3_con,DEG_cluster_4_con,DEG_cluster_5_con,DEG_cluster_6_con,DEG_cluster_7_con)
rownames(DEG_cluster_con)=DEG_cluster_con[,1]
Node_gene=DEG_cluster_con[NetGene,]
Node_gene=Node_gene[,-1]
Node_gene=na.omit(Node_gene)
Netsize=dim(Node_gene)[1] 

###### Correlation

PPI=PPI[rownames(Node_gene),rownames(Node_gene)]
sum(PPI)  #=820  edges
setwd("F:/parasite RNAseq Data Analysis/Result")  
write.table(PPI, file="PPI_code.xls",col.names = NA,sep = "\t") 

PPC= matrix(NA,nrow=Netsize,ncol=Netsize)  # Partial Pearson correlation coefficient 
PPC_p= matrix(NA,nrow=Netsize,ncol=Netsize)  # p value 

PC= matrix(NA,nrow=Netsize,ncol=Netsize)  # Partial Pearson correlation coefficient 
PC_p= matrix(NA,nrow=Netsize,ncol=Netsize)  # p value



for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    
      A1=pcor.test(as.numeric(Node_gene[i,]),as.numeric(Node_gene[j,]),as.numeric(as.matrix(Node_gene[-c(i,j),])))  # ,use = c("mat"), method = "pearson", na.rm =F 
      PPC[i,j]=A1$estimate
      PPC_p[i,j]=A1$p.value
    
    B=cor.test(as.numeric(Node_gene[i,]),as.numeric(Node_gene[j,]),method="pearson")
    PC[i,j]=B$estimate
    PC_p[i,j]=B$p.value
  }
}


PC[is.na(PC)]=0
PC_p[is.na(PC_p)]=1


PPC[is.na(PPC)]=0
PPC_p[is.na(PPC_p)]=1



#####  Save PPC 
setwd("F:/parasite RNAseq Data Analysis/Result")  
write.table(PPC, file="Partial correlation.xls",col.names = NA,sep = "\t") 
write.table(PPC_p, file="p value of Partial correlation.xls",col.names = NA,sep = "\t") 

### calculate the edge of the network

PPC=read.table("F:/parasite RNAseq Data Analysis/Result/Partial correlation.txt")
PPC_p=read.table("F:/parasite RNAseq Data Analysis/Result/p value of Partial correlation.txt")
Partial_PC=PPC
Partial_PC_p=PPC_p

##### Go to Matlab code Edge_CorrelationNetwork 

# edge=matrix(0,nrow=Netsize,ncol=Netsize)
# rownames(edge)=rownames(Node_gene); colnames(edge)=rownames(Node_gene); 
# 
# for (i in 1:Netsize){
#   for (j in 1:Netsize){
#     
#     if ((Partial_PC[i,j]>0.6) & (Partial_PC_p[i,j]<0.01) & (PPI[i,j]==1))
#       {edge[i,j]=1}
#     else 
#     {     
#       if ((Partial_PC[i,j]<-0.6) & (Partial_PC_p[i,j]<0.01) & (PPI[i,j]==1))
#     {edge[i,j]=-1}
#       else
#       {edge[i,j]=0}
#     }
#   }
# }
# 
# 
# sum(sum(edge!=0))  # 201

edge=read.csv("F:/parasite RNAseq Data Analysis/Result/edge_CorrelationNetwoek.csv",header=F)  # read the result from Matlab
rownames(edge)=rownames(Node_gene); colnames(edge)=rownames(Node_gene); 
edge["Serping1","Fn1"]
edge["Fn1","Serping1"]
edge["Serping1","Oasl2"]
edge["Mmp3","Oasl2"]  
  
### Select genes in the network
A=rownames(edge[rowSums(abs(edge))!=0, ])
B=colnames(edge[,colSums(abs(edge))!=0])
Net_Gene=union(A,B)

edge=edge[Net_Gene,Net_Gene]


dim(edge)  # 133,133


write.table(edge, file="edge matrix.txt",col.names = NA,sep = "\t") 


#### save gene expression data in the initial network
Node_gene=DEG_cluster_con[Net_Gene,]


write.table(Node_gene, file="Node.xls",col.names = NA,sep = "\t") 


