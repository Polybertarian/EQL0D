function [MAT,nubar,leak] = updateRates(MAT,RUN,RR)
    %MAT = UPDATERATES(MAT,RUN,RR) Updates reaction rates in case of PCC
    coeffs = interpCoeffs(RUN);
    for i=1:numel(MAT) %SYS.IDX.MAT.inFlux
        MAT(i).mCaptXS = 0.0*MAT(i).mCaptXS;
        MAT(i).mFissXS = 0.0*MAT(i).mFissXS;
        MAT(i).mN2nXS  = 0.0*MAT(i).mN2nXS;
        MAT(i).mN3nXS  = 0.0*MAT(i).mN3nXS;
        for j=1:3
            MAT(i).mCaptXS = MAT(i).mCaptXS + coeffs(j)*RR(j).devCapt(i)*RR(j).inMat.capt;
            MAT(i).mFissXS = MAT(i).mFissXS + coeffs(j)*RR(j).devFiss(i)*RR(j).inMat.fiss;
            MAT(i).mN2nXS  = MAT(i).mN2nXS  + coeffs(j)*RR(j).inMat.n2n;
            MAT(i).mN3nXS  = MAT(i).mN3nXS  + coeffs(j)*RR(j).inMat.n3n;
        end
    end
    nubar=coeffs*vertcat(RR.NU);
    leak=coeffs*vertcat(RR.LEAK);
end
