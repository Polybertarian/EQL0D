function [srcMat,dstMat] = transferNuclides(srcMat,dstMat,nucIdx,nucDens)
    % TRANSFERNUCLIDES(SRCMAT,DSTMAT,NUCIDX,NUCDENS) transfers nucDens nuclides identified by nucIdx from srcMat to dstMat
    if any(srcMat.N(nucIdx,end)<nucDens)&&all(nucDens)>0
        error('Error: Nuclide density of source material insufficient for transfer')
    elseif any(dstMat.N(nucIdx,end)<-nucDens)&&all(nucDens)<0
        error('Error: Nuclide density of destination material insufficient for transfer')
    else
        if ~isempty(srcMat)
            srcMat.N(nucIdx,end)=srcMat.N(nucIdx,end)-nucDens;
        end
        if ~isempty(dstMat)
            dstMat.N(nucIdx,end)=dstMat.N(nucIdx,end)+nucDens;
        end
    end
end
