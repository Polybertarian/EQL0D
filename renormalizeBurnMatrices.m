function SYS = renormalizeBurnMatrices(MAT,SYS)
%RENORMALIZERATES Takes burnup matrices in the system and renormalizes the
%flux-dependant rates by renormFactor (previously calculated) 

%%% Compare current rate to target
currentRate=0;
for i=SYS.IDX.fluxMat
    currentRate=currentRate+MAT(i).totReactionRate(SYS.RR.inMat.fiss);
end

fprintf(SYS.FID.log,'%s\n',['** RENORM ** Current flux: ' num2str([SYS.averageFlux],'%E')]);
SYS.renormFactor=SYS.targetFissionRate./currentRate;
SYS.averageFlux=SYS.averageFlux*SYS.renormFactor;
fprintf(SYS.FID.log,'%s\n',['** RENORM ** Current fission rate: ' num2str(currentRate,'%E') ...
    ', expected: ' num2str(SYS.targetFissionRate,'%E') ', factor: ' num2str(SYS.renormFactor)]);
if(SYS.renormFactor<0)
    error('Error: Renormalization factor negative! Something went wrong.')
end

%%% Change flux in materials
for j=SYS.IDX.fluxMat
    MAT(j).averageFlux=MAT(j).averageFlux*SYS.renormFactor;
    if(ismember(j,SYS.IDX.burnMat))
        SYS.MTX.burn{j}=SYS.MTX.burn{j}*SYS.renormFactor;
    end
end

%%% Change flux in non-included materials
for i=1:length(SYS.RR.notInMat)
   SYS.RR.notInMat(i).fiss=SYS.renormFactor*SYS.RR.notInMat(i).fiss; 
   SYS.RR.notInMat(i).capt=SYS.renormFactor*SYS.RR.notInMat(i).capt; 
   SYS.RR.notInMat(i).n2n=SYS.renormFactor*SYS.RR.notInMat(i).n2n; 
   SYS.RR.notInMat(i).n3n=SYS.renormFactor*SYS.RR.notInMat(i).n3n; 
end

currentRate=0;
for i=SYS.IDX.fluxMat
    currentRate=currentRate+MAT(i).totReactionRate(SYS.RR.inMat.fiss);
end

fprintf(SYS.FID.log,'%s\n',['** RENORM ** New fission rate: ' num2str(currentRate,'%E')]);

%%% Sum matrices
SYS.burnMatrix=[];
for i=SYS.IDX.contMat
    SYS.burnMatrix=blkdiag(SYS.burnMatrix,SYS.MTX.burn{j,i}+SYS.MTX.decay{j,i});
end
for i=1:length(SYS.MTX.rep(j,:))
    SYS.burnMatrix=SYS.burnMatrix+SYS.MTX.rep{j,i};
end

return
end

