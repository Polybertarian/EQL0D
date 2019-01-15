function [MAT,SYS] = burnCycle(MAT,OPT,REP,SYS)
    %[MAT,SYS] = burnCycle(MAT,OPT,REP,SYS) depletes the materials in the System

    if OPT.renormalize
        [MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat] = renormalizeSystem(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.tgtFissRate); % Renormalize burn matrices to new fission rate
        SYS = buildSystemMatrices(MAT,REP,SYS); % Build global matrix
    end

    for j=1:length(SYS.IDX.REP.matGroups)
        burnVec{j}=[]; % Prepare vector for solving
        for k=SYS.IDX.REP.matGroups{j}
            burnVec{j}=vertcat(burnVec{j},MAT(k).atDens(MAT(k).burnIdx));
        end
    end
    for j=1:length(SYS.IDX.REP.matGroups)
        if SYS.RUN.PCC.corrector || ~SYS.RUN.PCC.active
            tStep=SYS.RUN.tStep(end)/OPT.nSubSteps;
            burnVec2{j}=CRAMsolve(SYS.MTX.total{2,j},tStep,burnVec{j});
            for i=1:1:OPT.nSubSteps-1
                burnVec2{j}=CRAMsolve(SYS.MTX.total{2,j},tStep,burnVec2{j});
            end
        else
            tStep=SYS.RUN.tStep(end)*OPT.nSteps(SYS.RUN.ouCntr)/OPT.nSubSteps;
            burnVec2{j}=CRAMsolve(SYS.MTX.total{2,j},tStep,burnVec{j});
            for i=1:1:OPT.nSubSteps-1
                burnVec2{j}=CRAMsolve(SYS.MTX.total{2,j},tStep,burnVec2{j});
            end
        end
    end
    for j=1:length(SYS.IDX.REP.matGroups)
        for i=SYS.IDX.REP.matGroups{j} % Change material compositions
            MAT(i).N(MAT(i).burnIdx,end+1)=burnVec2{j}(SYS.MTX.burnMat{2,j}==i)*MAT(i).volume;
        end
    end

    if SYS.RUN.PCC.corrector
        prefix='C'; prefix2='CORR ';
        SYS.RUN.nowTime(end+1)=SYS.RUN.nowTime(end)+SYS.RUN.tStep(end)/(24.0*3600.0);
    else
        prefix='P'; prefix2='PRED ';
        if SYS.RUN.PCC.active
            SYS.RUN.nowTime(end+1)=SYS.RUN.nowTime(end)+SYS.RUN.tStep(end)*OPT.nSteps(SYS.RUN.ouCntr)/(24.0*3600.0); % in EFPD
        else
            SYS.RUN.nowTime(end+1)=SYS.RUN.nowTime(end)+SYS.RUN.tStep(end)/(24.0*3600.0); % in EFPD
        end
    end

    if strcmp(OPT.iterMode,'steps')||SYS.RUN.inCntr==1||floor(SYS.RUN.inCntr/10)==ceil(SYS.RUN.inCntr/10)&&~isempty(SYS.IDX.REP.batch)
        printStatus(MAT(SYS.IDX.MAT.burn),SYS.RUN,prefix2,'EoS (BB)');            % Add status to log
    end

    [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK); % Compute k-eff
    printK(SYS,'BB',prefix,'EQL0D');

    if OPT.renormalize
        [MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat] = renormalizeSystem(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.tgtFissRate);
    end

    if OPT.printSteps&&OPT.printStepsBB
        for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
            MAT(i).printMaterial(SYS.RUN.ouCntr,SYS.RUN.inCntr,SYS.RUN.nowTime,'BB'); % Print EOS composition
        end
    end

    SYS.prevN.EOC=[]; % Store compositions from previous loop
    for i=1:length(MAT)
        SYS.prevN.EOC=[SYS.prevN.EOC MAT(i).N(:,end)];
    end

    for k=1:numel(MAT) % Save compositions before batch
        MAT(k).N(:,end+1)=MAT(k).N(:,end);
    end
    SYS.RUN.nowTime(end+1)=SYS.RUN.nowTime(end); % Update time

    if OPT.redoxControl  % Adjust redox
        for i=SYS.IDX.redoxMat
            [MAT(i),dN] = MAT(i).redoxControl(SYS.IDX.redoxHalide{i},SYS.IDX.redoxNuc{i},OPT.REDOX.replaceMode);
            %name = MAT(i).nuclideName(SYS.IDX.redoxHalide{i});
            fprintf('%s\n',['** REDOX ** Excess of ' num2str(dN(1)*1E24,'%E') ' valence corrected.']);
        end
    end

    if ~isempty(SYS.IDX.REP.batch)  % Batchwise EoS processes
        if SYS.verboseMode
            fprintf('%s\n','** BATCH ** Performing batch processing steps...');
        end
        for r=SYS.IDX.REP.batch % Loop on batch processing streams
            if SYS.verboseMode
                fprintf('%s\n',['** BATCH ** Performing processing step ' REP(r).name '...']);
            end
            [MAT(REP(r).dstMatIdx),MAT(REP(r).srcMatIdx)] = REP(r).batchProcessing(MAT(REP(r).dstMatIdx),...
            MAT(REP(r).srcMatIdx),SYS.RUN.tStep(end));
        end
        if SYS.verboseMode
            fprintf('%s\n','** BATCH ** Batch processing steps finished!');
        end
    end
    [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)] = computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK); % Compute k-eff

    if OPT.reactControl  % Adjust reactivity
        [MAT,SYS] = reactivityControl(MAT,SYS,SYS.REA,SYS.IDX.REA);
    end

    %% After Batch processing
    if strcmp(OPT.iterMode,'steps')||SYS.RUN.inCntr==1||floor(SYS.RUN.inCntr/10)==ceil(SYS.RUN.inCntr/10)
        printStatus(MAT(SYS.IDX.MAT.burn),SYS.RUN,prefix2,'EoS (AB)');  % Add status to log
    end
    if SYS.RUN.inCntr<OPT.nSteps(SYS.RUN.ouCntr)
        [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK); % Compute k-eff
        printK(SYS,'AB',prefix,'EQL0D');
    end

    if OPT.renormalize
        [MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat] = renormalizeSystem(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.tgtFissRate);
    end

    if OPT.printSteps     % Print material composition to txt file
        for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
            MAT(i).printMaterial(SYS.RUN.ouCntr,SYS.RUN.inCntr,SYS.RUN.nowTime,'AB');
        end
    end
    for i=SYS.IDX.MAT.burn
        MAT(i).write(OPT.matWriteStyle); % Write compositions to serp files
    end
end
