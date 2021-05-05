function [A,B]=WNMFCDA(inputdata,params)

% Weighted Nonnegative Matrix Factorization based on  
%  Multi-Source Fusion Information for Predicting CircRNA-Disease Associations. 
%  2021/04/30

K=params.K; 
a=params.a;
p=params.p;
CD_mat=inputdata.CD_mat;
SC_mat=inputdata.SC_mat;
SD_mat=inputdata.SD_mat;


M = WKNN( CD_mat, SC_mat, SD_mat, K, a )
C_graph_mat = Graph( SC_mat , p );
D_graph_mat = Graph( SD_mat , p ); 
S_C = SC_mat.* C_graph_mat ;
S_D = SD_mat.* D_graph_mat ;  

clear K a p;
 
k=params.k;
lmada=params.lmada; 
beta=params.beta;
iterate=params.iterate; 
beta > 0  && lamda >0;
fprintf('k=%d  maxiter=%d  beta=%d  lamda=%d\n', k, iterate, beta, lamda ); 

[rows,cols] = size(M);
A=abs(rand(k,rows));        
B=abs(rand(k,cols));

D_c = diag(sum(S_C,2));
D_d = diag(sum(S_D,2));
L_c=D_c-S_C;
L_d=D_d-S_D;


fid = fopen( 'RunResult.txt','wt+');
for step=1:iterate
        U1=A.*((B*M'+lmada*A*S_C)./(B*B'*A+lmada*A*D_c+beta*A));
        V1=B.*((U1*M+lmada*B*S_D)./(U1*U1'*B+lmada*B*D_d+beta*B));
         
        ULU = sum(diag((U1*L_c)*U1'));
        VLV = sum(diag((V1*L_d)*V1'));
        obj = sum(sum((M-U1'*V1).^2))+beta*(sum(sum(U1.^2)) )+beta*(sum(sum(V1.^2)))+lmada*ULU+lmada*VLV; 
        
        error=max([max(sum((U1-A).^2)),max(sum((V1-B).^2))]);      
        
        fprintf(fid,'%s\n',[sprintf('step = \t'),int2str(step),...
            sprintf('\t obj = \t'),num2str(obj),...
		    sprintf('\t error = \t'), num2str(error)]);
        fprintf('step=%d  obj=%d  error=%d\n',step, obj, error);   
        if error< 10^(-4)
            fprintf('step=%d\n',step);
            break;
        end
        
        A=U1; 
        B=V1;
        
end
fclose(fid);

end

 
