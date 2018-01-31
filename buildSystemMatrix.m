function SYS = buildSystemMatrix(SYS)
%BUILDSYSTEMMATRIX Builds final burn-up matrix from all burn, decay and
%reprocessing matrices

%%% Remove old matrix
SYS.MTX.total(1)=[];

%%% Add new by summing burn, decay and rep matrices
SYS.MTX.total{2}=[];
for i=[SYS.IDX.burnMat SYS.IDX.contStrMat]
    SYS.MTX.total{2}=blkdiag(SYS.MTX.total{2},SYS.MTX.burn{2,i}+SYS.MTX.decay{2,i});
end
for i=SYS.IDX.contStr
    SYS.MTX.total{2}=SYS.MTX.total{2}+SYS.MTX.rep{2,i};
end

return
end

