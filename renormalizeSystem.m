function [MAT,SYS] = renormalizeSystem(MAT,SYS)
%RENORMALIZESYSTEM Takes burnup matrices in the system and renormalizes the
%flux-dependant rates by renormFactor (previously calculated)
global FID

coeffs = interpCoeffs(SYS);

%%% Compare current rate to target
currentRate=sum(SYS.RR(3).notInMat.fiss)*coeffs(3);
for j=1:2
    currentRate=currentRate+sum(SYS.RR(j).notInMat.fiss)*coeffs(j);
end
for i=SYS.IDX.MAT.inFlux
    currentRate=currentRate+sum(MAT(i).fissRate);
end

SYS.renormFactor=SYS.tgtFissRate./currentRate;
    fprintf('%s\n',['** RENORM ** Current fission rate: ' num2str(currentRate,'%E') ...
        ', expected: ' num2str(SYS.tgtFissRate,'%E') ', factor: ' num2str(SYS.renormFactor)]);
if SYS.renormFactor<0
    error('Error: Renormalization factor negative! Something went wrong.')
end
SYS.intFlux=SYS.intFlux*SYS.renormFactor;
for i=SYS.IDX.MAT.inFlux
    MAT(i).intFlux=MAT(i).intFlux*SYS.renormFactor;
end

%%% Change flux in non-included materials
for i=1:length(SYS.RR(2).notInMat)
    SYS.RR(2).notInMat(i).fiss=SYS.renormFactor*SYS.RR(2).notInMat(i).fiss;
    SYS.RR(2).notInMat(i).capt=SYS.renormFactor*SYS.RR(2).notInMat(i).capt;
    SYS.RR(2).notInMat(i).n2n=SYS.renormFactor*SYS.RR(2).notInMat(i).n2n;
    SYS.RR(2).notInMat(i).n3n=SYS.renormFactor*SYS.RR(2).notInMat(i).n3n;
end

currentRate=sum(SYS.RR(3).notInMat.fiss);
for i=SYS.IDX.MAT.inFlux
    currentRate=currentRate+sum(MAT(i).fissRate);
end
return
end
