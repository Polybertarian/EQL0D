function [MAT,SYS] = burnCycle(MAT,OPT,REP,SYS)
%[MAT,SYS] = burnCycle(MAT,OPT,REP,SYS) depletes the materials in the System

if OPT.renormalize
    [MAT,SYS] = renormalizeSystem(MAT,SYS); % renormalize burn matrices to new fission rate
    SYS = buildSystemMatrices(MAT,REP,SYS);  % build global matrix
end

parfor j=1:length(SYS.IDX.REP.matGroups)
    burnVec{j}=[]; %%% Prepare vector for solving
    for k=SYS.IDX.REP.matGroups{j}
        burnVec{j}=vertcat(burnVec{j},MAT(k).atDens(MAT(k).burnIdx));
    end
end
parfor j=1:length(SYS.IDX.REP.matGroups)
    if SYS.PCC.corrector || ~SYS.PCC.active
        tStep=SYS.tStep(end)/OPT.nSubSteps;
        burnVec2{j}=CRAMsolve(SYS.MTX.total{2,j},tStep,burnVec{j});
        for i=1:1:OPT.nSubSteps-1
            burnVec2{j}=CRAMsolve(SYS.MTX.total{2,j},tStep,burnVec2{j});
        end
    else
        tStep=SYS.tStep(end)*OPT.nSteps(SYS.ouCntr)/OPT.nSubSteps;
        burnVec2{j}=CRAMsolve(SYS.MTX.total{2,j},tStep,burnVec{j});
        for i=1:1:OPT.nSubSteps-1
            burnVec2{j}=CRAMsolve(SYS.MTX.total{2,j},tStep,burnVec2{j});
        end
    end
end
for j=1:length(SYS.IDX.REP.matGroups)
    for i=SYS.IDX.REP.matGroups{j} %%% Change material compositions
        MAT(i).N(MAT(i).burnIdx,end+1)=burnVec2{j}(SYS.MTX.burnMat{2,j}==i)*MAT(i).volume;
    end
end

if SYS.PCC.corrector
    prefix='C'; prefix2='CORRECTOR';
    SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0);
else
    prefix='P'; prefix2='PREDICTOR';
    if SYS.PCC.active
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)*OPT.nSteps(SYS.ouCntr)/(24.0*3600.0); %in EFPD
    else
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0); %in EFPD
    end
end

if strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10)&&~isempty(SYS.IDX.REP.batch)
    printStatus(MAT,SYS,prefix2,'EoS (BB)');            %%% Add status to log
end

[SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT,SYS); %%% Compute k-eff
printK(SYS,'BB',prefix,'EQL0D');

if OPT.renormalize
    [MAT,SYS] = renormalizeSystem(MAT,SYS);
end

if OPT.printSteps&&OPT.printStepsBB
    for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
        MAT(i).printMaterial(SYS.ouCntr,SYS.inCntr,SYS.nowTime,'BB'); %%% Print EOS composition
    end
end

SYS.prevN.EOC=[];
for i=1:length(MAT)
    SYS.prevN.EOC=[SYS.prevN.EOC MAT(i).N(:,end)]; %%% Store compositions from previous loop
end

%% Redox control
if OPT.redoxControl  %%% Adjust redox
    for i=SYS.IDX.redoxMat
        [MAT(i),dN] = MAT(i).redoxControl(SYS.IDX.redoxHalide{i},SYS.IDX.redoxNuc{i},OPT.REDOX.replaceMode);
        %name = MAT(i).nuclideName(SYS.IDX.redoxHalide{i});
        fprintf('%s\n',['** REDOX ** Excess of ' num2str(dN(1)*1E24,'%E') ' valence corrected.']);
    end
end

%% Batch processes
if ~isempty(SYS.IDX.REP.batch)  %%% Batchwise EoS processes
    if SYS.verboseMode
        fprintf('%s\n','** BATCH ** Performing batch processing steps...');
    end
    for r=SYS.IDX.REP.batch %%% Loop on batch processing streams
        if SYS.verboseMode
            fprintf('%s\n',['** BATCH ** Performing processing step ' REP(r).name '...']);
        end
        [MAT(REP(r).dstMatIdx),MAT(REP(r).srcMatIdx)] = REP(r).batchProcessing(MAT(REP(r).dstMatIdx),...
            MAT(REP(r).srcMatIdx),SYS.tStep(end));
    end
    if SYS.verboseMode
        fprintf('%s\n','** BATCH ** Batch processing steps finished!');
    end
end
[SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)] = computeK(MAT,SYS); %%% Compute k-eff

%% Reactivity control
if OPT.reactControl  %%% Adjust reactivity
    [MAT,SYS] = reactivityControl(MAT,SYS);
end

%% After Batch processing
if strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10)
    printStatus(MAT,SYS,prefix2,'EoS (AB)');  %%% Add status to log
end
if SYS.inCntr<OPT.nSteps(SYS.ouCntr)
    [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT,SYS); %%% Compute k-eff
    printK(SYS,'AB',prefix,'EQL0D');
end

if OPT.renormalize
    [MAT,SYS] = renormalizeSystem(MAT,SYS);
end

%%% Print material composition to file
if OPT.printSteps
    for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
        MAT(i).printMaterial(SYS.ouCntr,SYS.inCntr,SYS.nowTime,'AB');
    end
end
parfor i=SYS.IDX.MAT.burn
    MAT(i).write(OPT.matWriteStyle); %%% Write compositions
end
end
