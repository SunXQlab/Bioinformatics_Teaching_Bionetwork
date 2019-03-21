

##### read DEGs 

install.packages("ggm")
library(ggm)

source("pcor.R")

if(!require("ppcor")){
  install.packages("ppcor", repos='http://cran.us.r-project.org')
  library(ppcor)
}


setwd("F:/教学/生物信息学/生物信息学 上课 资料/")
DEG=read.table("DEGs.xls",fill= FALSE)


##### Initial PPI network
PPI_pair=read.table("PPI_interactions.txt",header = FALSE)
NetGene12=PPI_pair[,1:2]
NetGene=union(NetGene12[,1],NetGene12[,2])


Node_gene=DEG[NetGene,]
Node_gene=na.omit(Node_gene) 
Netsize=dim(Node_gene)[1] 
  
PPI=matrix(0,nrow=Netsize,ncol=Netsize)
for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    # for (k in 1:dim(NetGene12)[1])
    # {
    #   if (sum(c(NetGene[i],NetGene[j])==as.matrix(NetGene12)[k,])==2)  # 
    #     PPI[i,j]=PPI[i,j]+1
    #   else
    #     PPI[i,j]=PPI[i,j]+0  
    # }
    a=as.character(which(NetGene12[,1]==NetGene[i]))
    b=as.character(which(NetGene12[,2]==NetGene[j]))
    # a[is.null(a)]=0
    # b[is.null(b)]=0
    
    if (length(a)==0){a=0}
    if (length(b)==0){b=0}
    
    
    if (length(intersect(a,b))>0)
      PPI[i,j]=1
    else 
      PPI[i,j]=0
  }  
}


# for (i in 1:Netsize)
# {
#   for (j in 1:Netsize)
#   {
#     PPI[i,j]=max(PPI[i,j],PPI[j,i])
#   }
# }

rownames(PPI)=rownames(Node_gene)
colnames(PPI)=rownames(Node_gene)
sum(PPI)  #=830  edges




# 
# rownames(DEG)=DEG[,1]
# Node_gene=DEG[NetGene,]
# Node_gene=Node_gene[,-1]
# Node_gene=na.omit(Node_gene)
# Netsize=dim(Node_gene)[1] 

###### Correlation

PC= matrix(NA,nrow=Netsize,ncol=Netsize)  # Partial Pearson correlation coefficient 
PC_p= matrix(NA,nrow=Netsize,ncol=Netsize)  # p value



for (i in 1:Netsize)
{
  for (j in 1:Netsize)
  {
    
    B=cor.test(as.numeric(Node_gene[i,]),as.numeric(Node_gene[j,]),method="pearson")
    PC[i,j]=B$estimate
    PC_p[i,j]=B$p.value
  }
}


### calculate the edge of the network




edge=matrix(0,nrow=Netsize,ncol=Netsize)
rownames(edge)=rownames(Node_gene); colnames(edge)=rownames(Node_gene);

for (i in 1:Netsize){
  for (j in 1:Netsize){

    if ((abs(PC[i,j])>0.6) & (PC_p[i,j]<0.01) & (PPI[i,j]==1))
      {edge[i,j]=PC[i,j]}
    else
    {
      {edge[i,j]=0}
    }
  }
}

# 
# sum(sum(edge!=0))  # 201


### Reform edge data for cytoscape 

C=matrix(NA,nrow=sum(sum(edge!=0)),ncol=3)
row=1
for (i in 1:dim(edge)[1])
{
  for (j in 1:dim(edge)[1])
  {
    if (edge[i,j]!=0)
    {
      C[row,1]=rownames(edge)[i]
      C[row,3]=colnames(edge)[j]
      C[row,2]=edge[i,j]
      row=row+1
    }
  }
}

write.table(C, file="Correlation Network.txt",quote=F,row.names=F,col.names = F,sep = "\t")

