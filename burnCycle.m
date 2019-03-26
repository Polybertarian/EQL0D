function [MAT,SYS] = burnCycle(MAT,OPT,REP,SYS)
%[MAT,SYS] = burnCycle(MAT,OPT,REP,SYS) depletes the materials in the System

if OPT.renormalize
    [MAT,SYS] = renormalizeSystem(MAT,SYS); % Renormalize burn matrices to new fission rate
    SYS = buildSystemMatrices(MAT,REP,SYS); % Build global matrix
end

parfor j=1:length(SYS.IDX.REP.matGroups) % Prepare vector for solving
    burnVec{j}=[];
    for k=SYS.IDX.REP.matGroups{j}
        burnVec{j}=vertcat(burnVec{j},MAT(k).atDens(MAT(k).burnIdx));
    end
end

parfor j=1:length(SYS.IDX.REP.matGroups) % Solve system
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

for j=1:length(SYS.IDX.REP.matGroups) % Change material compositions
    for i=SYS.IDX.REP.matGroups{j}
        MAT(i).N(MAT(i).burnIdx,end+1)=burnVec2{j}(SYS.MTX.burnMat{2,j}==i)*MAT(i).volume;
    end
end

if SYS.PCC.corrector % Change time
    prefix='C'; prefix2='CORRECTOR';
    SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0); % in EFPD
else
    prefix='P'; prefix2='PREDICTOR';
    if SYS.PCC.active
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)*OPT.nSteps(SYS.ouCntr)/(24.0*3600.0); % in EFPD
    else
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0); % in EFPD
    end
end

if strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10)&&~isempty(SYS.IDX.REP.batch)
    printStatus(MAT,SYS,prefix2,'EoS (BB)');  % Add status to log
end

[SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT,SYS); % Compute k-eff
printK(SYS,'BB',prefix,'EQL0D');

if OPT.renormalize % Renormalize reaction rate to fission rate
    [MAT,SYS] = renormalizeSystem(MAT,SYS);
end

if OPT.printSteps&&OPT.printStepsBB % Print EOS composition
    for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
        MAT(i).printMaterial(SYS,'BB');
    end
end

SYS.prevN.EOC=[]; % Store compositions from previous loop
for i=1:length(MAT)
    SYS.prevN.EOC=[SYS.prevN.EOC MAT(i).N(:,end)];
end

if ~isempty(SYS.IDX.REP.batch)  % Batchwise end-of-step processes
    if SYS.verboseMode
        fprintf('%s\n','** BATCH ** Performing batch processing steps...');
    end
    for r=SYS.IDX.REP.batch % Loop on batch processing streams
        if SYS.verboseMode
            fprintf('%s\n',['** BATCH ** Performing processing step ' REP(r).name '...']);
        end
        [MAT(REP(r).dstMatIdx),MAT(REP(r).srcMatIdx)] = REP(r).batchProcessing(MAT(REP(r).dstMatIdx),MAT(REP(r).srcMatIdx),SYS.tStep(end));
    end
    if SYS.verboseMode
        fprintf('%s\n','** BATCH ** Batch processing steps finished!');
    end
end

if OPT.redoxControl  % Adjust fluorine/chlorine concentration
    for i=SYS.IDX.redoxMat
        [MAT(i),dN] = MAT(i).redoxControl(SYS.IDX.redoxHalide{i},SYS.IDX.redoxNuc{i},OPT.REDOX.replaceMode);
        name = MAT(i).nuclideName(SYS.IDX.redoxHalide{i});
        fprintf('%s\n',['** REDOX ** Excess of ' num2str(dN(1)*1E24,'%E') ' valence corrected.']);
    end
end

[SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)] = computeK(MAT,SYS); % Compute k-eff

if OPT.reactControl  % Reactivity control
    [MAT,SYS] = reactivityControl(MAT,OPT.REA,SYS);
end

if strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10) % Add status to log
    printStatus(MAT,SYS,prefix2,'EoS (AB)');
end

if SYS.inCntr<OPT.nSteps(SYS.ouCntr) % Compute k-eff
    [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT,SYS);
    printK(SYS,'AB',prefix,'EQL0D');
end

if OPT.renormalize
    [MAT,SYS] = renormalizeSystem(MAT,SYS);
end

if OPT.printSteps % Print material composition to file
    for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
        MAT(i).printMaterial(SYS,'AB');
    end
end
parfor i=SYS.IDX.MAT.burn % Write compositions to Serpent files
    MAT(i).write(OPT.matWriteStyle);
end

end
