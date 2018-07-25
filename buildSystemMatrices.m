function SYS = buildSystemMatrices(MAT,REP,SYS)
%BUILDSYSTEMMATRICES Builds final burn-up matrix from all burn, decay and
%reprocessing matrices


SYS.MTX.total(1,:)=[];
for i=1:length(SYS.IDX.REP.matGroups)
  SYS.MTX.total{2,i}=[];
  SYS.MTX.burnMat{2,i}=[];
  streams=[];
  for j=SYS.IDX.REP.matGroups{i}
    SYS.MTX.total{2,i}=blkdiag(SYS.MTX.total{2,i},MAT(j).normBurnMtx+MAT(j).decMtx);
    SYS.MTX.burnMat{2,i}=vertcat(SYS.MTX.burnMat{2,i},j*ones(length(MAT(j).burnZAI),1));
    streams=unique([streams MAT(j).streams.cont]);
  end
  for k=streams
    globalRepMtx=0.0*SYS.MTX.total{2,i};
    idxSrc=SYS.MTX.burnMat{2,i}==REP(k).srcMatIdx;
    idxDst=SYS.MTX.burnMat{2,i}==REP(k).dstMatIdx;
    if(REP(k).isKeep)
      globalRepMtx(idxSrc,idxSrc)=createKeepMatrix(REP(k),MAT(REP(k).srcMatIdx),MAT(REP(k).dstMatIdx));
    else
      globalRepMtx(idxSrc,idxSrc)=REP(k).repMtx;
    end
    globalRepMtx(idxSrc,idxDst)=-(MAT(REP(k).srcMatIdx).volume/MAT(REP(k).dstMatIdx).volume)*globalRepMtx(idxSrc,idxSrc);
    SYS.MTX.total{2,i}=SYS.MTX.total{2,i}+globalRepMtx;
  end
end

return
end

