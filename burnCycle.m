function [MAT,SYS] = burnCycle(MAT,OPT,REP,SYS)
%[MAT,SYS] = burnCycle(MAT,OPT,REP,SYS) depletes the materials in the
%System
%% Burning
if(SYS.PCC.corrector)
    prefix='C';
    prefix2='CORRECTOR';
    if(SYS.inCntr==1)
        SYS.nowTime(end+1)=SYS.nowTime(end)-SYS.tStep(end)*OPT.nSteps(SYS.ouCntr)/(24.0*3600.0);
        for i=[SYS.IDX.burnMat SYS.IDX.contStrMat]
            MAT(i).N(:,end+1)=SYS.prevN.BOC(:,i);
        end
    end
    MAT=updateRates(MAT,SYS);
    SYS.MTX.interp=((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.MTX.total{1}+...
        (SYS.inCntr-1)*SYS.MTX.total{2})/SYS.PCC.nSteps; %%% Interpolate and solve
    SYS.burnVec=CRAMsolve(SYS.MTX.interp,SYS.tStep(end),SYS.burnVec);
    SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0);
    for i=[SYS.IDX.burnMat SYS.IDX.contStrMat] %%% Change material compositions
        MAT(i).N(SYS.IDX.burnZAI{2,i},end+1)=SYS.burnVec(vertcat(SYS.IDX.matZAI{2,:})==i)*MAT(i).volume;
    end
else
    prefix='P';
    prefix2='PREDICTOR';
    SYS.burnVec=[]; %%% Prepare vector for solving
    for i=[SYS.IDX.burnMat SYS.IDX.contStrMat]
        SYS.burnVec=vertcat(SYS.burnVec,MAT(i).atDens(SYS.IDX.burnZAI{2,i}));
    end
    if(SYS.PCC.active)
        SYS.burnVec2=CRAMsolve(SYS.MTX.total{2},SYS.tStep(end)*OPT.nSteps(SYS.ouCntr),SYS.burnVec);
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)*OPT.nSteps(SYS.ouCntr)/(24.0*3600.0); %in EFPD
    else
        SYS.burnVec2=CRAMsolve(SYS.MTX.total{2},SYS.tStep(end),SYS.burnVec);
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0); %in EFPD
    end
    for i=[SYS.IDX.burnMat SYS.IDX.contStrMat] %%% Change material compositions
        MAT(i).N(SYS.IDX.burnZAI{2,i},end+1)=SYS.burnVec2(vertcat(SYS.IDX.matZAI{2,:})==i)*MAT(i).volume;
    end
end

%% Before Batch Processing

if(strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10))
    printStatus(MAT,SYS,prefix2,'EoS (BB)');            %%% Add status to log
end

SYS=computeK(MAT,SYS); %%% Compute k-eff
printK(SYS,'BB',prefix,'EQL0D');

if(OPT.printSteps&&OPT.printStepsBB)
    for i=[SYS.IDX.burnMat SYS.IDX.strMat]
        MAT(i).printMaterial(SYS,'BB'); %%% Print EOS composition
    end
end

SYS.prevN.EOC=[];
for i=1:length(MAT)
    SYS.prevN.EOC=[SYS.prevN.EOC MAT(i).N(:,end)]; %%% Store compositions from previous loop
end

if(~OPT.PCC||SYS.PCC.corrector)
    if(OPT.redoxControl) %%% Adjust redox
        for i=SYS.IDX.redoxMat
            [MAT(i),dN] = MAT(i).redoxControl(SYS.IDX.redoxHalide{i},SYS.IDX.redoxNuc{i},OPT.REDOX.replaceMode);
            name=MAT.nuclideName(SYS.IDX.redoxHalide{i});
            fprintf(SYS.FID.log,'%s\n',['** REDOX ** Excess of ' num2str(-dN*1E24,'%E') ' ' [name{:}]  ' atoms corrected.']);
        end
    end
    if(~isempty(SYS.IDX.batchStr)) %%% Batchwise EoS processes
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
                MAT(SYS.IDX.srcMat(r)).N(:,end+1)=MAT(SYS.IDX.srcMat(r)).N(:,end);
                MAT(SYS.IDX.srcMat(r)).N(SYS.IDX.srcPos{r},end)=MAT(SYS.IDX.srcMat(r)).N(SYS.IDX.srcPos{r},end)-share.*mDefect;
            end
            if(SYS.IDX.dstMat(r)~=0)
                MAT(SYS.IDX.dstMat(r)).N(:,end+1)=MAT(SYS.IDX.dstMat(r)).N(:,end);
                MAT(SYS.IDX.dstMat(r)).N(SYS.IDX.dstPos{r},end)=MAT(SYS.IDX.dstMat(r)).N(SYS.IDX.dstPos{r},end)+share.*mDefect;
            end
            clearvars mDefect share
        end
        if(SYS.verboseMode)
            fprintf(SYS.FID.log,'%s\n','** BATCH ** Batch processing steps finished!');
        end
    end
    SYS=computeK(MAT,SYS); %%% Compute k-eff
    if(OPT.reactControl) %%% Adjust reactivity
        [MAT,SYS] = reactivityControl(MAT,OPT.REA,SYS);
    end
end
%% After Batch processing
if(strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10))
    printStatus(MAT,SYS,prefix2,'EoS (AB)');  %%% Add status to log
end

if(SYS.inCntr<OPT.nSteps(SYS.ouCntr))
    SYS=computeK(MAT,SYS); %%% Compute k-eff
    printK(SYS,'AB',prefix,'EQL0D');
end

%%% Print material composition to file
if(OPT.printSteps)
    for i=[SYS.IDX.burnMat SYS.IDX.strMat]
        MAT(i).printMaterial(SYS,'AB');
    end
end

for i=SYS.IDX.burnMat
    MAT(i).write(OPT.matWriteStyle); %%% Write compositions
end

end

