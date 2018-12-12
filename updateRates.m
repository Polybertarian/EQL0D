function MAT = updateRates(MAT,SYS)
%MAT = UPDATERATES(MAT,SYS) Updates reaction rates in case of PCC

coeffs = interpCoeffs(SYS);

for i=SYS.IDX.MAT.inFlux
    MAT(i).mCaptXS = coeffs(3)*SYS.RR(3).devCapt(i)*SYS.RR(3).inMat.capt...
        +coeffs(2)*SYS.RR(2).devCapt(i)*SYS.RR(2).inMat.capt...
        +coeffs(1)*SYS.RR(1).devCapt(i)*SYS.RR(1).inMat.capt;
    MAT(i).mFissXS = coeffs(3)*SYS.RR(3).devFiss(i)*SYS.RR(3).inMat.fiss...
        +coeffs(2)*SYS.RR(2).devFiss(i)*SYS.RR(2).inMat.fiss...
        +coeffs(1)*SYS.RR(1).devFiss(i)*SYS.RR(1).inMat.fiss;
    MAT(i).mN2nXS  = coeffs(3)*SYS.RR(3).inMat.n2n +coeffs(2)*SYS.RR(2).inMat.n2n...
        +coeffs(1)*SYS.RR(1).inMat.n2n;
    MAT(i).mN3nXS  = coeffs(3)*SYS.RR(3).inMat.n3n +coeffs(2)*SYS.RR(2).inMat.n3n...
        +coeffs(1)*SYS.RR(1).inMat.n3n;
end

return
end

