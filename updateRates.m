function MAT = updateRates(MAT,SYS)
%UPDATERATES Updates reaction rates in case of PCC

if(SYS.PCC.active)
    for i=SYS.IDX.fluxMat
        MAT(i).mCaptXS=((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.RR.inMat{1}.capt+...
            (SYS.inCntr)*SYS.RR.inMat{2}.capt)/SYS.PCC.nSteps;
        MAT(i).mFissXS=((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.RR.inMat{1}.fiss+...
            (SYS.inCntr)*SYS.RR.inMat{2}.fiss)/SYS.PCC.nSteps;
        MAT(i).mN2nXS =((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.RR.inMat{1}.n2n+... 
            (SYS.inCntr)*SYS.RR.inMat{2}.n2n)/SYS.PCC.nSteps;
        MAT(i).mN3nXS =((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.RR.inMat{1}.n3n+...
            (SYS.inCntr)*SYS.RR.inMat{2}.n3n)/SYS.PCC.nSteps;
    end
else
    for i=SYS.IDX.fluxMat
        MAT(i).mCaptXS=SYS.RR.inMat{2}.capt;
        MAT(i).mFissXS=SYS.RR.inMat{2}.fiss;
        MAT(i).mN2nXS =SYS.RR.inMat{2}.n2n;
        MAT(i).mN3nXS =SYS.RR.inMat{2}.n3n;
    end
end
return
end

