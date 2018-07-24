function SYS = buildSystemMatrix(MAT,SYS)
%BUILDSYSTEMMATRIX Builds final burn-up matrix from all burn, decay and
%reprocessing matrices

SYS.MTX.burn={};
for i=1:length(SYS.IDX.REP.matGroups)
  SYS.MTX.burn{i}=[];
  SYS.MTX.burnMat{i}=[];
  streams=[];
  for j=SYS.IDX.REP.matGroups{i}
    SYS.MTX.burn{i}=blkdiag(SYS.MTX.burn{i},MAT(j).normBurnMtx+MAT(j).decMtx);
    SYS.MTX.burnMat{i}=vertcat(SYS.MTX.burnMat{i},j*ones(length(MAT(j).burnZAI),1));
    streams=unique([streams MAT(j).streams.cont]);
  end
  for k=streams
    globalRepMtx=0.0*SYS.MTX.burn{i};
    idxSrc=SYS.MTX.burnMat{i}==REP(k).srcMatIdx;
    idxDst=SYS.MTX.burnMat{i}==REP(k).dstMatIdx;
    if(REP(k).isKeep)
      
    else
      globalRepMtx(idxSrc,idxSrc)=REP(k).repMtx;
      globalRepMtx(idxSrc,idxDst)=-(MAT(REP(k).srcMatIdx).volume/MAT(REP(k).dstMatIdx).volume)*REP(k).repMtx;
    end
    SYS.MTX.burn{i}=SYS.MTX.burn{i}+globalRepMtx;
  end
end

%%% Remove old matrix
SYS.MTX.total(1)=[];

%%% Add new by summing burn, decay and rep matrices
SYS.MTX.total{2}=[];
for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
  SYS.MTX.total{2}=blkdiag(SYS.MTX.total{2},SYS.MTX.burn{2,i}+SYS.MTX.decay{2,i});
end
for i=SYS.IDX.REP.cont
  SYS.MTX.total{2}=SYS.MTX.total{2}+SYS.MTX.rep{2,i};
end

return
end

