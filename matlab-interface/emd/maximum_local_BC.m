% function [Index,Value]=maximum_local_BC(Matrix)
% traitement des données aux bords : création d'une périodicité dans les
% deux directions (ce qui sortz à droite re-rentre à gauche, ce qui sort en 
% haut re-rentre en bas et vice versa)
% i.e. si matrice (2:mx+1)*(2:nx+1), [1,:]=[mx,:] [mx+2,:]=[3,:]
% [:,1]=[:,nx] [:,nx+2]=[:,3]

function [Index,Value]=maximum_local_BC(Matrix)

[m,n]=size(Matrix);
Index=[0 0];
Value=0;
Matrix(2:m+1,2:n+1)=Matrix;
Matrix(1,2:n+1)=Matrix(m,2:n+1);
Matrix(m+2,2:n+1)=Matrix(3,2:n+1);
Matrix(:,1)=Matrix(:,n);
Matrix(:,n+2)=Matrix(:,3);

for i=2:m+1
    for j=2:n+1
        greater=true(1);
        if     Matrix(i,j)<=Matrix(i-1,j-1),  greater=false(1);
        elseif Matrix(i,j)<=Matrix(i-1,j),    greater=false(1);
        elseif Matrix(i,j)<=Matrix(i-1,j+1),  greater=false(1);
        elseif Matrix(i,j)<=Matrix(i,j-1),    greater=false(1);
        elseif Matrix(i,j)<=Matrix(i,j+1),    greater=false(1);
        elseif Matrix(i,j)<=Matrix(i+1,j-1),  greater=false(1);
        elseif Matrix(i,j)<=Matrix(i+1,j),    greater=false(1);
        elseif Matrix(i,j)<=Matrix(i+1,j+1),  greater=false(1);
        end
        if greater==true(1) Index=[Index;i-1 j-1]; Value=[Value;Matrix(i,j)];end
    end    
end

% delete the first line of Index and Value
dim=size(Index);
Index=Index(2:dim(1),:);
Value=Value(2:length(Value));