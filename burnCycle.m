function [MAT,SYS] = burnCycle(MAT,OPT,REP,SYS)
%[MAT,SYS] = burnCycle(MAT,OPT,REP,SYS) depletes the materials in the
%System
%% Burning
SYS.burnVec=[]; %%% Prepare vector for solving
for i=[SYS.IDX.MAT.burn SYS.IDX.REP.contMat]
    SYS.burnVec=vertcat(SYS.burnVec,MAT(i).atDens(SYS.IDX.burnZAI{2,i}));
end
if(SYS.PCC.corrector)
    prefix='C';
    prefix2='CORRECTOR';
    [MAT,SYS.MTX.interp]=updateRates(MAT,SYS);
    SYS.burnVec2=CRAMsolve(SYS.MTX.interp,SYS.tStep(end),SYS.burnVec);
    SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0);
else
    prefix='P';
    prefix2='PREDICTOR';
    if(SYS.PCC.active)
        SYS.burnVec2=CRAMsolve(SYS.MTX.total{2},SYS.tStep(end)*OPT.nSteps(SYS.ouCntr),SYS.burnVec);
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)*OPT.nSteps(SYS.ouCntr)/(24.0*3600.0); %in EFPD
    else
        SYS.burnVec2=CRAMsolve(SYS.MTX.total{2},SYS.tStep(end),SYS.burnVec);
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0); %in EFPD
    end
end
for i=[SYS.IDX.MAT.burn SYS.IDX.REP.contMat] %%% Change material compositions
	MAT(i).N(SYS.IDX.burnZAI{2,i},end+1)=SYS.burnVec2(vertcat(SYS.IDX.matZAI{2,:})==i)*MAT(i).volume;
end

%% Before Batch Processing
if(strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10))
    printStatus(MAT,SYS,prefix2,'EoS (BB)');            %%% Add status to log
end

[SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT,SYS); %%% Compute k-eff
printK(SYS,'BB',prefix,'EQL0D');

if(OPT.printSteps&&OPT.printStepsBB)
    for i=[SYS.IDX.MAT.burn SYS.IDX.strMat]
        MAT(i).printMaterial(SYS,'BB'); %%% Print EOS composition
    end
end

SYS.prevN.EOC=[];
for i=1:length(MAT)
    SYS.prevN.EOC=[SYS.prevN.EOC MAT(i).N(:,end)]; %%% Store compositions from previous loop
end

%if(~OPT.PCC||SYS.PCC.corrector)
    if(OPT.redoxControl) %%% Adjust redox
        for i=SYS.IDX.redoxMat
            [MAT(i),dN] = MAT(i).redoxControl(SYS.IDX.redoxHalide{i},SYS.IDX.redoxNuc{i},OPT.REDOX.replaceMode);
            name=MAT(i).nuclideName(SYS.IDX.redoxHalide{i});
            fprintf(SYS.FID.log,'%s\n',['** REDOX ** Excess of ' num2str(dN(1)*1E24,'%E') ' valence corrected.']);
        end
    end
    if(~isempty(SYS.IDX.REP.batch)) %%% Batchwise EoS processes
        %%% Performs batch-wise processing steps
        if(SYS.verboseMode)
            fprintf(SYS.FID.log,'%s\n','** BATCH ** Performing batch processing steps...');
        end
        %%% Loop on batch processing streams
        for r=SYS.IDX.REP.batch
            if(SYS.verboseMode)
                fprintf(SYS.FID.log,'%s\n',['** BATCH ** Performing processing step ' REP(r).name '...']);
            end
            [MAT(REP(r).dstMatIdx),MAT(REP(r).srcMatIdx)] = ...
                REP(r).batchProcessing(MAT(REP(r).dstMatIdx),MAT(REP(r).srcMatIdx),SYS.tStep(end));
        end
        if(SYS.verboseMode)
            fprintf(SYS.FID.log,'%s\n','** BATCH ** Batch processing steps finished!');
        end
    end
    [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT,SYS); %%% Compute k-eff
    if(OPT.reactControl) %%% Adjust reactivity
        [MAT,SYS] = reactivityControl(MAT,OPT.REA,SYS);
    end
%end
%% After Batch processing
if(strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10))
    printStatus(MAT,SYS,prefix2,'EoS (AB)');  %%% Add status to log
end

if(SYS.inCntr<OPT.nSteps(SYS.ouCntr))
    [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT,SYS); %%% Compute k-eff
    printK(SYS,'AB',prefix,'EQL0D');
end

%%% Print material composition to file
if(OPT.printSteps)
    for i=[SYS.IDX.MAT.burn SYS.IDX.strMat]
        MAT(i).printMaterial(SYS,'AB');
    end
end

for i=SYS.IDX.MAT.burn
    MAT(i).write(OPT.matWriteStyle); %%% Write compositions
end

end

