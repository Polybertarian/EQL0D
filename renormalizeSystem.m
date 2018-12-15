function [fluxMAT,notInMatRR] = renormalizeSystem(fluxMAT,notInMatRR,targetFissRate)
%RENORMALIZESYSTEM Takes burnup matrices in the system and renormalizes the
%flux-dependant rates by renormFactor (previously calculated)

currentRate=sum(notInMatRR.fiss);
for i=1:numel(fluxMAT)%SYS.IDX.MAT.inFlux
    currentRate=currentRate+sum(fluxMAT(i).fissRate);
end

renormFactor=targetFissRate./currentRate; % Compare current rate to target
fprintf('%s\n',['** RENORM ** Current fission rate: ' num2str(currentRate,'%E') ...
    ', expected: ' num2str(targetFissRate,'%E') ', factor: ' num2str(renormFactor)]);
if renormFactor<0
    error('Error: Renormalization factor negative! Something went wrong.')
end

for i=1:numel(fluxMAT)%SYS.IDX.MAT.inFlux
    fluxMAT(i).intFlux=fluxMAT(i).intFlux*renormFactor;
end

notInMatRR.fiss=renormFactor*notInMatRR.fiss; % Change flux in non-included materials
notInMatRR.capt=renormFactor*notInMatRR.capt;
notInMatRR.n2n=renormFactor*notInMatRR.n2n;
notInMatRR.n3n=renormFactor*notInMatRR.n3n;
end
