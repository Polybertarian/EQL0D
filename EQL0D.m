function EQL0D(SYS)
    %EQL0D(SYS) is the main function of the EQL0D procedure, containing the
    %outer and inner loops
    global FID
    try
        [MAT,OPT,REP,SYS] = initialize(SYS); % Initialize or restart from .mat file
        while ~SYS.stopOuter % Outer Loop / Cycles
            SYS.RUN.ouCntr=SYS.RUN.ouCntr+1; SYS.RUN.inCntr=0;
            if strcmp(OPT.iterMode,'equilibrium')
                SYS.RUN.tStep(end+1)=OPT.cycleLength(SYS.RUN.ouCntr)*24.0*3600.0;
            elseif strcmp(OPT.iterMode,'steps')
                SYS.RUN.tStep(end+1)=OPT.cycleLength(SYS.RUN.ouCntr)*24.0*3600.0/OPT.nSteps(SYS.RUN.ouCntr);
            end
            SYS.stopInner=false; SYS.RUN.PCC.active=OPT.PCC; SYS.RUN.PCC.corrector=false;SYS.oldFIMA=[MAT.FIMA];
            SYS.oldN=[];
            for i=1:length(MAT)
                SYS.oldN=[SYS.oldN MAT(i).N(:,end)]; %%% Store compositions from previous loop
            end
            if ~SYS.debugMode
                modifyInput(SYS.Casename,'dep',OPT.cycleLength(SYS.RUN.ouCntr)); %%% Adapt depletion time/burnup in Serpent
                runSerpent(OPT,SYS); % run serpent
            end
            if ~SYS.debugMode||SYS.RUN.ouCntr==1
                [MAT,SYS] = loadSerpentData(MAT,SYS); %%% Read Serpent outputs
                if ~SYS.debugMode
                    if SYS.RUN.ouCntr==1
                        SYS.RUN.ouCntr=0;
                    end
                    for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                        if OPT.printSteps
                            MAT(i).printMaterial(SYS,'AB'); % write material composition to txt file
                        elseif OPT.printCycles&&~OPT.printSteps
                            MAT(i).printMaterial(SYS,'EoC'); % write material composition to txt file
                        end
                    end
                    if SYS.RUN.ouCntr==0
                        SYS.RUN.ouCntr=1;
                    end
                    saveFiles({MAT(SYS.IDX.MAT.burn).name},OPT.filesToKeep,SYS);  %%% Move files to folder
                    printK(SYS,'AB','C','Serpent'); % print keff and kinf to file
                end
            end
            [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK); %%% Compute k-eff
            printK(SYS,'AB','C','EQL0D'); % print keff and kinf to file
            if OPT.renormalize
                [MAT,SYS] = renormalizeSystem(MAT,SYS);  % renormalize burn matrices to new fission rate
            end
            SYS = buildSystemMatrices(MAT,REP,SYS); % build global matrix
            save([SYS.Casename '.mat']);  % save to .mat file
            while ~SYS.stopInner % Inner loop
                SYS.RUN.inCntr=SYS.RUN.inCntr+1; SYS.prevFIMA=[MAT(SYS.IDX.MAT.burn).FIMA];
                SYS.prevN.BOC=[];
                for i=1:length(MAT)
                    SYS.prevN.BOC=[SYS.prevN.BOC MAT(i).N(:,end)]; % Store compositions from previous loop
                end

                [MAT,SYS]=burnCycle(MAT,OPT,REP,SYS);  % deplete materials and perform batch operations

                if OPT.PCC&&~SYS.debugMode %Corrector step
                    SYS.RUN.PCC.corrector=true; SYS.RUN.PCC.nSteps=OPT.nSteps(SYS.RUN.ouCntr);
                    runSerpent(OPT,SYS); % Run Serpent/Read outputs
                    [MAT,SYS] = loadSerpentData(MAT,SYS);
                    SYS.RUN.nowTime(end+1)=SYS.RUN.nowTime(end)-SYS.RUN.tStep(end)*OPT.nSteps(SYS.RUN.ouCntr)/(24.0*3600.0);
                    for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                        MAT(i).N(:,end+1)=SYS.prevN.BOC(:,i);
                    end
                    if ~isempty(SYS.IDX.REP.cont)
                        SYS = createKeepMatrices(MAT,REP,SYS);
                    end
                    SYS = buildSystemMatrices(MAT,REP,SYS);
                    if ~SYS.debugMode
                        saveFiles({MAT(SYS.IDX.MAT.burn).name},OPT.filesToKeep,SYS);
                        printK(SYS,'AB','P','Serpent');
                    end
                    [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK);
                    printK(SYS,'AB','P','EQL0D');

                    while SYS.RUN.inCntr<=OPT.nSteps(SYS.RUN.ouCntr) %%% Burn system
                        [MAT,SYS]=burnCycle(MAT,OPT,REP,SYS);
                        SYS.RUN.inCntr=SYS.RUN.inCntr+1;
                    end
                    SYS.RUN.PCC.corrector=false; SYS.stopInner=true;
                else
                    SYS=testConvergence(MAT,OPT,SYS,'inner');
                end
                %save([SYS.Casename '.mat']);
            end
            if OPT.printCycles&&~OPT.printSteps
                for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                    MAT(i).printMaterial(SYS,'EoC'); %%% Print compositions at EoC
                end
            end
            SYS=testConvergence(MAT,OPT,SYS,'outer'); %%% stop outer loop?
        end

        %% Last step / Equilibrium results
        if ~SYS.debugMode
            runSerpent(OPT,SYS);
        end
        [MAT,SYS] = loadSerpentData(MAT,SYS);
        switch OPT.iterMode
        case 'equilibrium'
            for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                MAT(i).printMaterial(SYS,'EQL_AB');
                mat=MAT(i); mat.N(:,end+1)=SYS.prevN.EOC(:,i);
                mat.printMaterial(SYS,'EQL_BB');
            end
            neutronBalance(MAT,SYS)
        case 'steps'
            if OPT.printSteps
                for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                    MAT(i).printMaterial(SYS,'EoC');
                end
            end
            if ~SYS.debugMode
                SYS.RUN.ouCntr=SYS.RUN.ouCntr+1;
                saveFiles({MAT(SYS.IDX.MAT.burn).name},OPT.filesToKeep,SYS);
            end
        end
        save([SYS.Casename '.mat']);
        fprintf('%s\n','**** EQL0D **** Procedure finished.');
        for file=fields(FID)
            fclose(FID.(file{1}));
        end
        if OPT.writeMail
            unix(['echo "...in ' pwd ' !" | mail -s "EQL0D calculation finished!" $LOGNAME'])
        end
        return
    catch exception %%% error handling
        fprintf('%s\n',getReport(exception,'extended'));
        save([SYS.Casename '_err.mat'])
        [~,~]=unix(['echo "...in ' pwd ' !" | mail -s "EQL0D crashed!" $LOGNAME']);
        exit(1)
    end
end
