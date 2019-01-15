function [srcMat,dstMat,nucDensRepl]=replaceNuclides(srcMat,dstMat,nucIdx,nucDens,replIdx,replFrac)
    % REPLACENUCLIDES(SRCMAT,DSTMAT,NUCIDX,NUCDENS,REPLIDX) transfers nucDens nuclides identified by nucIdx from srcMat to dstMat
    nucDensRepl=replFrac.*sum(nucDens.*dstMat.atomicMass(nucIdx))./dstMat.avMass(replIdx);
    if any(dstMat.N(replIdx,end)<abs(nucDensRepl))
        error('Error: Nuclide density of destination material insufficient for transfer')
    else
        [srcMat,dstMat] = transferNuclides(srcMat,dstMat,nucIdx,nucDens); %take nucDens from src to dst
        [dstMat,srcMat] = transferNuclides(dstMat,srcMat,replIdx,nucDensRepl); %take nucdens from dst to src
    end
end
