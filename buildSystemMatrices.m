function SYS = buildSystemMatrices(MAT,REP,SYS)
%BUILDSYSTEMMATRICES Builds final burn-up matrix from all burn, decay and
%reprocessing matrices

SYS.MTX.total={};
for i=1:length(SYS.IDX.REP.matGroups)
  SYS.MTX.total{i}=[];
  SYS.MTX.burnMat{i}=[];
  streams=[];
  for j=SYS.IDX.REP.matGroups{i}
    SYS.MTX.total{i}=blkdiag(SYS.MTX.total{i},MAT(j).normBurnMtx+MAT(j).decMtx);
    SYS.MTX.burnMat{i}=vertcat(SYS.MTX.burnMat{i},j*ones(length(MAT(j).burnZAI),1));
    streams=unique([streams MAT(j).streams.cont]);
  end
  for k=streams
    globalRepMtx=0.0*SYS.MTX.total{i};
    idxSrc=SYS.MTX.burnMat{i}==REP(k).srcMatIdx;
    idxDst=SYS.MTX.burnMat{i}==REP(k).dstMatIdx;
    if(REP(k).isKeep)
      globalRepMtx(idxSrc,idxSrc)=createKeepMatrix(REP(k),MAT(REP(k).srcMatIdx),MAT(REP(k).dstMatIdx));
    else
      globalRepMtx(idxSrc,idxSrc)=REP(k).repMtx;
    end
    globalRepMtx(idxSrc,idxDst)=-(MAT(REP(k).srcMatIdx).volume/MAT(REP(k).dstMatIdx).volume)*globalRepMtx(idxSrc,idxSrc);
    SYS.MTX.total{i}=SYS.MTX.total{i}+globalRepMtx;
  end
end

return
end

