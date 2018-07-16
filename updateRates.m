function [MAT, MTX] = updateRates(MAT,SYS)
%MAT = UPDATERATES(MAT,SYS) Updates reaction rates in case of PCC
if(SYS.PCC.corrector)
    coeffs=(SYS.inCntr-1)/SYS.PCC.nSteps;
    coeffs=[1-coeffs coeffs];
    mCaptXS=coeffs(1)*SYS.RR.inMat{1}.capt+coeffs(2)*SYS.RR.inMat{2}.capt;
    mFissXS=coeffs(1)*SYS.RR.inMat{1}.fiss+coeffs(2)*SYS.RR.inMat{2}.fiss;
    mN2nXS =coeffs(1)*SYS.RR.inMat{1}.n2n+coeffs(2)*SYS.RR.inMat{2}.n2n;
    mN3nXS =coeffs(1)*SYS.RR.inMat{1}.n3n+coeffs(2)*SYS.RR.inMat{2}.n3n;
    for i=SYS.IDX.MAT.inFlux
        MAT(i).mCaptXS=mCaptXS;
        MAT(i).mFissXS=mFissXS;
        MAT(i).mN2nXS=mN2nXS;
        MAT(i).mN3nXS=mN3nXS;
    end
    MTX=coeffs(1)*SYS.MTX.total{1}+coeffs(2)*SYS.MTX.total{2}; %%% Interpolate and solve 
else
    for i=SYS.IDX.MAT.inFlux
        MAT(i).mCaptXS=SYS.RR.inMat{2}.capt;
        MAT(i).mFissXS=SYS.RR.inMat{2}.fiss;
        MAT(i).mN2nXS =SYS.RR.inMat{2}.n2n;
        MAT(i).mN3nXS =SYS.RR.inMat{2}.n3n;
    end
    MTX=[];
end
return
end

