function [keff,kinf] = computeK(MAT,SYS)
%[keff,kinf] = COMPUTEK(MAT,SYS) compute k-eff/inf

coeffs = interpCoeffs(SYS);

nubar=coeffs(3)*SYS.RR(3).NU+coeffs(2)*SYS.RR(2).NU+coeffs(1)*SYS.RR(1).NU;
leak=coeffs(1)*SYS.RR(1).LEAK+coeffs(2)*SYS.RR(2).LEAK+coeffs(3)*SYS.RR(3).LEAK;

fiss=sum([SYS.RR(3).notInMat.fiss]); % usually 0 if all materials correctly input
prod=0;
capt=sum([SYS.RR(3).notInMat.capt]);
n2n =sum([SYS.RR(3).notInMat.n2n]);
n3n =sum([SYS.RR(3).notInMat.n3n]);

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
