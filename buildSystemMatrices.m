function SYS = buildSystemMatrices(MAT,REP,SYS)
%BUILDSYSTEMMATRICES Builds final burn-up matrix from all burn, decay and
%reprocessing matrices

coeffs = interpCoeffs(SYS);
SYS.MTX.total(1,:)=[];
for i=1:length(SYS.IDX.REP.matGroups)
    SYS.MTX.total{2,i}=[];
    SYS.MTX.burnMat{2,i}=[];
    streams=[];

    for j=SYS.IDX.REP.matGroups{i}
        SYS.MTX.total{2,i}=blkdiag(SYS.MTX.total{2,i},MAT(j).normBurnMtx(coeffs)+MAT(j).decMtx);
        SYS.MTX.burnMat{2,i}=vertcat(SYS.MTX.burnMat{2,i},j*ones(length(MAT(j).burnZAI),1));
        streams=unique([streams MAT(j).streams.cont]);
    end
    if ~isempty(streams)
        for k=streams
            globalRepMtx=0.0*SYS.MTX.total{2,i};
            if ~isempty(REP(k).srcMatIdx)
                idxSrc=SYS.MTX.burnMat{2,i}==REP(k).srcMatIdx;
            else
                idxSrc=[];
            end
            if ~isempty(REP(k).dstMatIdx)
                idxDst=SYS.MTX.burnMat{2,i}==REP(k).dstMatIdx;
            else
                idxDst=[];
            end
            if REP(k).isKeep
                globalRepMtx=createKeepMatrix(REP(k),MAT(REP(k).srcMatIdx),MAT(REP(k).dstMatIdx),...
                    SYS.MTX.total{2,i},SYS.MTX.burnMat{2,i});
            else
                globalRepMtx(idxSrc,idxSrc)=REP(k).repMtx;
            end
            if ~isempty(REP(k).dstMatIdx)&&~isempty(REP(k).srcMatIdx)
                if REP(k).isKeep
                else
                    globalRepMtx(idxDst,idxSrc)=-(MAT(REP(k).srcMatIdx).volume/MAT(REP(k).dstMatIdx).volume)*...
                        globalRepMtx(idxSrc,idxSrc);
                end
            end
            SYS.MTX.total{2,i}=SYS.MTX.total{2,i}+globalRepMtx;
        end
    end
end

return
end
