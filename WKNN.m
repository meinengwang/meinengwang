function [CD_mat_new] = WKNN( CD_mat, SC_mat, SD_mat, K, a )

% Weighted K nearest neighbor profiles for circRNAs and diseases

[rows,cols]=size(CD_mat);
M_c=zeros(rows,cols);  
M_d=zeros(rows,cols);  

knn_network_c = KNN( SC_mat, K );  %for circRNA
for i = 1 : rows   
         w=zeros(1,K);
        [sort_c,idx_c]=sort(knn_network_c(i,:),2,'descend'); 
        sum_c=sum(sort_c(1,1:K));   
        for j = 1 : K
            w(1,j)=a^(j-1)*sort_c(1,j); 
            M_c(i,:) =  M_c(i,:)+ w(1,j)* CD_mat(idx_c(1,j),:); 
        end                      
            M_c(i,:)=M_c(i,:)/sum_c;              
end

knn_network_d = KNN( SD_mat , K );  %for disease
for i = 1 : cols   
        w=zeros(1,K);
        [sort_d,idx_d]=sort(knn_network_d(i,:),2,'descend');
        sum_d=sum(sort_d(1,1:K));
        for j = 1 : K
            w(1,j)=a^(j-1)*sort_d(1,j);
            M_d(:,i) =  M_d(:,i)+ w(1,j)* CD_mat(:,idx_d(1,j)); 
        end                      
            M_d(:,i)=M_d(:,i)/sum_d;               
end

a1=1;
a2=1;
M_cd=(M_c*a1+M_d*a2)/(a1+a2);  

 for i = 1 : rows
     for j = 1 : cols
         CD_mat_new(i,j)=max(CD_mat(i,j),M_cd(i,j));
     end    
 end

end

function [ knn_network ] = KNN( network , k )
    [rows, cols] = size( network );
    network= network-diag(diag(network)); 
    knn_network = zeros(rows, cols);
    [sort_network,idx]=sort(network,2,'descend');
    for i = 1 : rows
        knn_network(i,idx(i,1:k))=sort_network(i,1:k);
    end
end


