function [keff,kinf] = computeK(fluxMAT,RUN,RR)
%[keff,kinf] = COMPUTEK(fluxMAT,RUN,RR) compute k-eff/inf

coeffs = interpCoeffs(RUN);

nubar=coeffs(3)*RR(3).NU+coeffs(2)*RR(2).NU+coeffs(1)*RR(1).NU;
leak=coeffs(1)*RR(1).LEAK+coeffs(2)*RR(2).LEAK+coeffs(3)*RR(3).LEAK;

prod=0;
fiss=sum([RR(3).notInMat.fiss]); % Usually 0 if all materials correctly input
capt=sum([RR(3).notInMat.capt]); n2n =sum([RR(3).notInMat.n2n]);
n3n =sum([RR(3).notInMat.n3n]);

for i=numel(fluxMAT) %SYS.IDX.MAT.inFlux
    fiss=fiss+sum(fluxMAT(i).fissRate); prod=prod+nubar(i)*sum(fluxMAT(i).fissRate);
    capt=capt+sum(fluxMAT(i).captRate); n2n =n2n+sum(fluxMAT(i).n2nRate);
    n3n =n3n+sum(fluxMAT(i).n3nRate);
end

kinf=(prod+2.0*n2n+3.0*n3n)/(fiss+capt+n2n+n3n); keff=kinf*leak;

return
end
