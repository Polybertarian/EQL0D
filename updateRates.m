function MAT = updateRates(MAT,SYS)
%MAT = UPDATERATES(MAT,SYS) Updates reaction rates in case of PCC

coeffs = interpCoeffs(SYS);

for i=SYS.IDX.MAT.inFlux
    MAT(i).mCaptXS=coeffs(3)*SYS.RR.devCapt(3,i)*SYS.RR.inMat{3}.capt...
        +coeffs(2)*SYS.RR.devCapt(2,i)*SYS.RR.inMat{2}.capt...
        +coeffs(1)*SYS.RR.devCapt(1,i)*SYS.RR.inMat{1}.capt;
    MAT(i).mFissXS=coeffs(3)*SYS.RR.devFiss(3,i)*SYS.RR.inMat{3}.fiss...
        +coeffs(2)*SYS.RR.devFiss(2,i)*SYS.RR.inMat{2}.fiss...
        +coeffs(1)*SYS.RR.devFiss(1,i)*SYS.RR.inMat{1}.fiss;
    MAT(i).mN2nXS =coeffs(3)*SYS.RR.inMat{3}.n2n +coeffs(2)*SYS.RR.inMat{2}.n2n...
        +coeffs(1)*SYS.RR.inMat{1}.n2n;
    MAT(i).mN3nXS =coeffs(3)*SYS.RR.inMat{3}.n3n +coeffs(2)*SYS.RR.inMat{2}.n3n...
        +coeffs(1)*SYS.RR.inMat{1}.n3n;
end

return

end
