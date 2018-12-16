function [srcMat,dstMat] = transferNuclides(srcMat,dstMat,nucIdx,nucDens)
% TRANSFERNUCLIDES(SRCMAT,DSTMAT,NUCIDX,NUCDENS) transfers nucDens nuclides identified by nucIdx from srcMat to dstMat
srcMat.N(nucIdx)=srcMat.N(nucIdx)-nucDens;
dstMat.N(nucIdx)=dstMat.N(nucIdx)+nucDens;
end
