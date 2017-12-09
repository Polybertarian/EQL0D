function volume=findVolume(Material,Rep,SYS,targetK)
fprintf(SYS.FID.log,'%s\n','Finding volume of feed to be added for criticality...');

vector=((SYS.RR.NU(end))*[SYS.RR.inMat.fiss]+2*[SYS.RR.inMat.n2n]+3*[SYS.RR.inMat.n3n]...
    -targetK*([SYS.RR.inMat.fiss]+[SYS.RR.inMat.n2n]+[SYS.RR.inMat.n3n]+[SYS.RR.inMat.capt]));

vector2=(SYS.RR.NU(end))*sum([SYS.RR.notInMat.fiss])+2*sum([SYS.RR.notInMat.n2n])+3*sum([SYS.RR.notInMat.n3n])...
    -targetK*(sum([SYS.RR.notInMat.fiss])+sum([SYS.RR.notInMat.n2n])+sum([SYS.RR.notInMat.n3n])+sum([SYS.RR.notInMat.capt]));

volume=-Material.volume*(vector2+vector*Material.atDens)/(vector2+vector(SYS.IDX.dstPos{SYS.IDX.addVolStr})*Rep.share);

return
end


