function [keff,kinf] = computeK(MAT,SYS)
%[keff,kinf] = COMPUTEK(MAT,SYS) compute k-inf and k-eff using reaction rates from Serpent and compositions
%from the MAT vector


coeffs = interpCoeffs(SYS);

nubar=coeffs(3)*SYS.RR.NU(3,:)+coeffs(2)*SYS.RR.NU(2,:)+coeffs(1)*SYS.RR.NU(1,:);
leak=coeffs(1)*SYS.RR.LEAK(1)+coeffs(2)*1*SYS.RR.LEAK(2)+coeffs(3)*1*SYS.RR.LEAK(3);

fiss=sum([SYS.RR.notInMat{end}.fiss]); %%% usually 0 if all materials correctly input
prod=0;
capt=sum([SYS.RR.notInMat{end}.capt]);
n2n =sum([SYS.RR.notInMat{end}.n2n]);
n3n =sum([SYS.RR.notInMat{end}.n3n]);

for i=SYS.IDX.MAT.inFlux
    fiss=fiss+sum(MAT(i).fissRate);
    prod=prod+nubar(i)*sum(MAT(i).fissRate);
    capt=capt+sum(MAT(i).captRate);
    n2n =n2n+sum(MAT(i).n2nRate);
    n3n =n3n+sum(MAT(i).n3nRate);
end

kinf=(prod+2.0*n2n+3.0*n3n)/(fiss+capt+n2n+n3n);
keff=kinf*leak;

return
end
