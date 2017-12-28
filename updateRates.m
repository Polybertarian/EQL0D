function MAT = updateRates(MAT,SYS)
%MAT = UPDATERATES(MAT,SYS) Updates reaction rates in case of PCC

if(SYS.PCC.corrector)
    mCaptXS=((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.RR.inMat{1}.capt+(SYS.inCntr)*SYS.RR.inMat{2}.capt)/SYS.PCC.nSteps;
    mFissXS=((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.RR.inMat{1}.fiss+(SYS.inCntr)*SYS.RR.inMat{2}.fiss)/SYS.PCC.nSteps;
    mN2nXS =((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.RR.inMat{1}.n2n+(SYS.inCntr)*SYS.RR.inMat{2}.n2n)/SYS.PCC.nSteps;
    mN3nXS =((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.RR.inMat{1}.n3n+(SYS.inCntr)*SYS.RR.inMat{2}.n3n)/SYS.PCC.nSteps;
    for i=SYS.IDX.fluxMat
        MAT(i).mCaptXS=mCaptXS;
        MAT(i).mFissXS=mFissXS;
        MAT(i).mN2nXS=mN2nXS;
        MAT(i).mN3nXS=mN3nXS;
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

