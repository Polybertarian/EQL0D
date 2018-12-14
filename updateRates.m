function MAT = updateRates(MAT,RUN,RR)
%MAT = UPDATERATES(MAT,SYS) Updates reaction rates in case of PCC

coeffs = interpCoeffs(RUN);

for i=1:numel(MAT) %SYS.IDX.MAT.inFlux
    MAT(i).mCaptXS = coeffs(3)*RR(3).devCapt(i)*RR(3).inMat.capt...
        +coeffs(2)*RR(2).devCapt(i)*RR(2).inMat.capt...
        +coeffs(1)*RR(1).devCapt(i)*RR(1).inMat.capt;
    MAT(i).mFissXS = coeffs(3)*RR(3).devFiss(i)*RR(3).inMat.fiss...
        +coeffs(2)*RR(2).devFiss(i)*RR(2).inMat.fiss...
        +coeffs(1)*RR(1).devFiss(i)*RR(1).inMat.fiss;
    MAT(i).mN2nXS  = coeffs(3)*RR(3).inMat.n2n +coeffs(2)*RR(2).inMat.n2n...
        +coeffs(1)*RR(1).inMat.n2n;
    MAT(i).mN3nXS  = coeffs(3)*RR(3).inMat.n3n +coeffs(2)*RR(2).inMat.n3n...
        +coeffs(1)*RR(1).inMat.n3n;
end

return
end
