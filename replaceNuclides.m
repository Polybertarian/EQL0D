function [srcMat,dstMat,nucDensRepl]=replaceNuclides(srcMat,dstMat,nucIdx,nucDens,replIdx,replFrac)
    % REPLACENUCLIDES(SRCMAT,DSTMAT,NUCIDX,NUCDENS,REPLIDX) transfers nucDens nuclides identified by nucIdx from srcMat to dstMat
    if any(dstMat.avMass(replIdx)==0)
        error('zero average mass detected before transfer!')
    end
    nucDensRepl=replFrac.*sum(nucDens.*dstMat.atomicMass(nucIdx))./dstMat.avMass(replIdx);
    if any(dstMat.N(replIdx,end)<nucDensRepl)
        error('Error: Nuclide density of destination material insufficient for transfer')
    else
        [srcMat,dstMat] = transferNuclides(srcMat,dstMat,nucIdx,nucDens);
        if any(dstMat.avMass(replIdx)==0)
            error('zero average mass detected after transfer!')
        end

        if any(isnan(nucDensRepl))
            error('NaN value detected!')
        end
        [dstMat,srcMat] = transferNuclides(dstMat,srcMat,replIdx,nucDensRepl);
    end
end
