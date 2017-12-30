function MAT = batchProcessing(MAT,REP,SYS)
%%% Performs batch-wise processing steps
if(SYS.verboseMode)
    fprintf(SYS.FID.log,'%s\n','** BATCH ** Performing batch processing steps...');
end
%%% Loop on batch processing streams
for r=SYS.IDX.batchStr
    if(SYS.verboseMode)
        fprintf(SYS.FID.log,'%s\n',['** BATCH ** Performing processing step ' REP(r).name '...']);
    end
    %%% 1)determine isotopics, 2)determine total quantity, 3) make change
    if(SYS.IDX.srcMat(r)==0)
        share=REP(r).share./MAT(SYS.IDX.dstMat(r)).atomicMass(SYS.IDX.dstPos{r})/1.0E24;
    else
        if(ismember(REP(r).mode,{'keepAFPM','keepAM','keepTotM'}))
            share=MAT(SYS.IDX.srcMat(r)).mFrac(SYS.IDX.srcPos{r});
            share=share./MAT(SYS.IDX.srcMat(r)).atomicMass(SYS.IDX.srcPos{r})/1.0E24;
        elseif(ismember(REP(r).mode,{'keepAFPA','keepAA','keepTotA'}))
            share=MAT(SYS.IDX.srcMat(r)).aFrac(SYS.IDX.srcPos{r});
        else
            share=MAT(SYS.IDX.srcMat(r)).aFrac(SYS.IDX.srcPos{r});
        end
    end
    switch REP(r).mode
        case 'keepTotM' %refill total mass
            mDefect=MAT(SYS.IDX.dstMat(r)).initTotMass-MAT(SYS.IDX.dstMat(r)).totMass;
        case 'keepTotA' %refill total nuclides
            mDefect=MAT(SYS.IDX.dstMat(r)).initTotN-MAT(SYS.IDX.dstMat(r)).totN;
        case 'keepAFPM' %refill up to initial actinide mass - current FP mass
            mDefect=MAT(SYS.IDX.dstMat(r)).initTotActMass-MAT(SYS.IDX.dstMat(r)).totActMass-MAT(SYS.IDX.dstMat(r)).totFPMass;
        case 'keepAFPA' %refill up to initial actinide mass - current FP mass
            mDefect=MAT(SYS.IDX.dstMat(r)).initTotActN-MAT(SYS.IDX.dstMat(r)).totActN-MAT(SYS.IDX.dstMat(r)).totFPN;
        case 'keepAM'   %refill up to initial actinide mass
            mDefect=MAT(SYS.IDX.dstMat(r)).initTotActMass-MAT(SYS.IDX.dstMat(r)).totActMass;
        case 'keepAA'   %refill up to initial actinide nuclides
            mDefect=MAT(SYS.IDX.dstMat(r)).initTotActN-MAT(SYS.IDX.dstMat(r)).totActN;
        case 'remove'
            mDefect=REP(r).rate*SYS.tStep(end)*MAT(SYS.IDX.srcMat(r)).N(SYS.IDX.srcPos{r},end); %not a mass
            share=REP(r).share;
    end
    if(SYS.IDX.srcMat(r)~=0)
        MAT(SYS.IDX.srcMat(r)).N(SYS.IDX.srcPos{r},end+1)=MAT(SYS.IDX.srcMat(r)).N(SYS.IDX.srcPos{r},end)-share.*mDefect;
    end
    if(SYS.IDX.dstMat(r)~=0)
        MAT(SYS.IDX.dstMat(r)).N(SYS.IDX.dstPos{r},end+1)=MAT(SYS.IDX.dstMat(r)).N(SYS.IDX.dstPos{r},end)+share.*mDefect;
    end
    clearvars mDefect share
end
if(SYS.verboseMode)
    fprintf(SYS.FID.log,'%s\n','** BATCH ** Batch processing steps finished!');
end

end
