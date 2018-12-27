function [srcMat,dstMat] = transferNuclides(srcMat,dstMat,nucIdx,nucDens)
    % TRANSFERNUCLIDES(SRCMAT,DSTMAT,NUCIDX,NUCDENS) transfers nucDens nuclides identified by nucIdx from srcMat to dstMat
    srcMat.N(nucIdx,end)=srcMat.N(nucIdx,end)-nucDens;
    if ~isempty(dstMat)
        dstMat.N(nucIdx,end)=dstMat.N(nucIdx,end)+nucDens;
    end
end
