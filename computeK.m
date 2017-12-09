function SYS = computeK(MAT,SYS)
%COMPUTEK compute k-inf and k-eff using reaction rates from Serpent and compositions
%from the MAT vector

fiss=sum([SYS.RR.notInMat{2}.fiss]); %%% usually 0 if all materials correctly input
capt=sum([SYS.RR.notInMat{2}.capt]);
n2n =sum([SYS.RR.notInMat{2}.n2n]);
n3n =sum([SYS.RR.notInMat{2}.n3n]);

for i=SYS.IDX.fluxMat
    fiss=fiss+sum(MAT(i).fissRate);
    capt=capt+sum(MAT(i).captRate);
    n2n =n2n+sum(MAT(i).n2nRate);
    n3n =n3n+sum(MAT(i).n3nRate);
end

if(SYS.PCC.active)
    kinf=(((SYS.PCC.nSteps-SYS.inCntr)*SYS.RR.NU{1}+(SYS.inCntr)*SYS.RR.NU{2})/SYS.PCC.nSteps*fiss+2.0*n2n+3.0*n3n)/(fiss+capt+n2n+n3n);
    keff=kinf*((SYS.PCC.nSteps-SYS.inCntr)*1/SYS.RR.LEAK{1}+(SYS.inCntr)*1/SYS.RR.LEAK{2})/SYS.PCC.nSteps;   
else
    kinf=(SYS.RR.NU{2}*fiss+2.0*n2n+3.0*n3n)/(fiss+capt+n2n+n3n);
    keff=kinf/SYS.RR.LEAK{2};
end

SYS.keff(end+1)=keff;
SYS.kinf(end+1)=kinf;
return
end
