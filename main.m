function main(SYS)
    %main(SYS) is the main function of the EQL0D procedure, containing the
    %outer and inner loops

    try
        [MAT,OPT,REP,SYS] = initialize(SYS); % Initialize or restart from .mat file
        while ~SYS.stopOuter % Outer Loop / Cycles
            SYS.ouCntr=SYS.ouCntr+1; SYS.inCntr=0;
            if strcmp(OPT.iterMode,'equilibrium')
                SYS.tStep(end+1)=OPT.cycleLength(SYS.ouCntr)*24.0*3600.0;
            elseif strcmp(OPT.iterMode,'steps')
                SYS.tStep(end+1)=OPT.cycleLength(SYS.ouCntr)*24.0*3600.0/OPT.nSteps(SYS.ouCntr);
            end
            SYS.stopInner=false; SYS.PCC.active=OPT.PCC; SYS.PCC.corrector=false;SYS.oldFIMA=[MAT.FIMA];
            SYS.oldN=[];
            for i=1:length(MAT)
                SYS.oldN=[SYS.oldN MAT(i).N(:,end)]; % Store compositions from previous loop
            end
            if ~SYS.debugMode
                modifyInput(SYS.Casename,'dep',OPT.cycleLength(SYS.ouCntr)); % Adapt depletion time/burnup in Serpent
                runSerpent(OPT.serpentPath,SYS.Casename,SYS.nCores); % Run serpent
            end
            if ~SYS.debugMode||SYS.ouCntr==1
                [MAT,SYS] = loadSerpentData(MAT,SYS); % Read Serpent outputs
                if ~SYS.debugMode
                    if SYS.ouCntr==1
                        SYS.ouCntr=0;
                    end
                    for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                        if OPT.printSteps
                            MAT(i).printMaterial(SYS,'AB'); % Write material composition to txt file
                        elseif OPT.printCycles&&~OPT.printSteps
                            MAT(i).printMaterial(SYS,'EoC'); % Write material composition to txt file
                        end
                    end
                    if SYS.ouCntr==0
                        SYS.ouCntr=1;
                    end
                    saveFiles({MAT(SYS.IDX.MAT.burn).name},OPT.keepFiles,SYS.Casename,SYS.ouCntr,~SYS.PCC.corrector&SYS.PCC.active); % Move files to folder
                    printK(SYS,'AB','C','Serpent'); % Print keff and kinf to file
                end
            end
            [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT,SYS); % Compute k-eff
            printK(SYS,'AB','C','EQL0D'); % Print keff and kinf to file
            if OPT.renormalize
                [MAT,SYS] = renormalizeSystem(MAT,SYS);  % Renormalize burn matrices to new fission rate
            end
            SYS = buildSystemMatrices(MAT,REP,SYS); % Build global matrix
            save([SYS.Casename '.mat']);  % Save to .mat file
            while ~SYS.stopInner  % Inner loop
                SYS.inCntr=SYS.inCntr+1; SYS.prevFIMA=[MAT(SYS.IDX.MAT.burn).FIMA];
                SYS.prevN.BOC=[];
                for i=1:length(MAT)
                    SYS.prevN.BOC=[SYS.prevN.BOC MAT(i).N(:,end)]; % Store compositions from previous loop
                end

                [MAT,SYS]=burnCycle(MAT,OPT,REP,SYS);  % Deplete materials and perform batch operations

                if OPT.PCC&&~SYS.debugMode  %Corrector step
                    SYS.PCC.corrector=true; SYS.PCC.nSteps=OPT.nSteps(SYS.ouCntr);
                    runSerpent(OPT.serpentPath,SYS.Casename,SYS.nCores); % Run Serpent/Read outputs
                    [MAT,SYS] = loadSerpentData(MAT,SYS);
                    SYS.nowTime(end+1)=SYS.nowTime(end)-SYS.tStep(end)*OPT.nSteps(SYS.ouCntr)/(24.0*3600.0);
                    for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                        MAT(i).N(:,end+1)=SYS.prevN.BOC(:,i);
                    end
                    SYS = buildSystemMatrices(MAT,REP,SYS);
                    saveFiles({MAT(SYS.IDX.MAT.burn).name},OPT.keepFiles,SYS.Casename,SYS.ouCntr,~SYS.PCC.corrector&SYS.PCC.active);
                    printK(SYS,'AB','P','Serpent');
                    [SYS.KEFF.EQL0D(end+1),SYS.KINF.EQL0D(end+1)]=computeK(MAT,SYS);
                    printK(SYS,'AB','P','EQL0D');

                    while SYS.inCntr<=OPT.nSteps(SYS.ouCntr) % Burn system
                        [MAT,SYS]=burnCycle(MAT,OPT,REP,SYS);
                        SYS.inCntr=SYS.inCntr+1;
                    end
                    SYS.PCC.corrector=false; SYS.stopInner=true;
                else
                    SYS=testConvergence(MAT,OPT,SYS,'inner');
                end
                %save([SYS.Casename '.mat']);
            end
            if OPT.printCycles&&~OPT.printSteps
                for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                    MAT(i).printMaterial(SYS,'EoC'); % Print compositions at EoC
                end
            end
            SYS=testConvergence(MAT,OPT,SYS,'outer'); % Stop outer loop?
        end

        if ~SYS.debugMode % Last step / Equilibrium results
            runSerpent(OPT.serpentPath,SYS.Casename,SYS.nCores);
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
            if(OPT.printSteps)
                for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                    MAT(i).printMaterial(SYS,'EoC');
                end
            end
            SYS.ouCntr=SYS.ouCntr+1;
            saveFiles({MAT(SYS.IDX.MAT.burn).name},OPT.keepFiles,SYS.Casename,SYS.ouCntr,~SYS.PCC.corrector&SYS.PCC.active);
        end
        save([SYS.Casename '.mat']);
        fprintf('%s\n','**** EQL0D **** Procedure finished.');
        for file=fields(SYS.FID)
            fclose(SYS.FID.(file{1}));
        end
        if OPT.writeMail
            unix(['echo "...in ' pwd ' !" | mail -s "EQL0D calculation finished!" $LOGNAME'])
        end
    catch exception %%% error handling
        fprintf('%s\n',getReport(exception,'extended'));
        save([SYS.Casename '.mat'])
        [~,~]=unix(['echo "...in ' pwd ' !" | mail -s "EQL0D crashed!" $LOGNAME']);
        exit(1)
    end
end
