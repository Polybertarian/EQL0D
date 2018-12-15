function [keff,kinf] = computeK(fluxMAT,notInMatRR,nubar,leak)
%[keff,kinf] = COMPUTEK(fluxMAT,RUN,RR) compute k-eff/inf

production=0;
fiss=sum([notInMatRR.fiss]); % Usually 0 if all materials correctly input
capt=sum([notInMatRR.capt]);
n2n =sum([notInMatRR.n2n]);
n3n =sum([notInMatRR.n3n]);

for i=1:numel(fluxMAT) %SYS.IDX.MAT.inFlux
    fiss=fiss+sum(fluxMAT(i).fissRate);
    production=production+nubar(i)*sum(fluxMAT(i).fissRate);
    capt=capt+sum(fluxMAT(i).captRate);
    n2n = n2n+sum(fluxMAT(i).n2nRate);
    n3n = n3n+sum(fluxMAT(i).n3nRate);
end

kinf=(production+2.0*n2n+3.0*n3n)/(fiss+capt+n2n+n3n);
keff=kinf*leak;
end
