function [A2,B2,c] = resizeMatrix(A,a,B,b)
%RESIZEMATRIX Resize matrix A,B and index vectors a,b to same size

%%% resize index vector
if(iscolumn(a)&&iscolumn(b))
    c=unique(sort(vertcat(a,b)));
elseif(~iscolumn(a)&&~iscolumn(b))
    c=unique(sort(horzcat(a,b)));
end
m=numel(c);

%%% get indexes
idxA2=ismember(c,a);
idxB2=ismember(c,b);

%%% resize matrix
if(issparse(A)&&issparse(B))
    A2=spalloc(m,m,numel(nonzeros(A)));
    B2=spalloc(m,m,numel(nonzeros(B)));
else
    A2=zeros(m,m);
    B2=zeros(m,m);
end
A2(idxA2,idxA2)=A;
B2(idxB2,idxB2)=B;

return
end

