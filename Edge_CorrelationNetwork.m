%%%  calculate the edge of the network
PPI=xlsread('F:\parasite RNAseq Data Analysis\Result\PPI_code.xls');

PPC=xlsread('F:\parasite RNAseq Data Analysis\Result\Partial correlation.xls');
PPC_p=xlsread('F:\parasite RNAseq Data Analysis/Result\p value of Partial correlation.xls');
Partial_PC=PPC(:,2:167);
Partial_PC_p=PPC_p(:,2:167);

edge=zeros(166,166);
Netsize=166;
% rownames(edge)=rownames(Node_gene); colnames(edge)=rownames(Node_gene); 

for i=1:Netsize
  for j =1:Netsize
    if ((Partial_PC(i,j)>0.9) && (Partial_PC_p(i,j)<0.05) && (PPI(i,j)==1))
      edge(i,j)=1;
    elseif ((Partial_PC(i,j)<-0.9) && (Partial_PC_p(i,j)<0.05) && (PPI(i,j)==1))
    edge(i,j)=-1;
    else
      edge(i,j)=0;
    end
  end
end

sum(sum(edge~=0))  % 460
xlswrite('F:\parasite RNAseq Data Analysis\Result\edge_CorrelationNetwoek.xls', edge);  

edge(60,25)  % to test the interaction between Fn1 and Serping1 that should be 1
edge(25,60)
edge(60,41)
edge(41,60)
edge(60,79)