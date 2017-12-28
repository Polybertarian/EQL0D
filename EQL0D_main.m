function [] = EQL0D_main(SYS)
%EQL0D_MAIN is the main function of the EQL0D procedure, containing the
%outer and inner loops

try
    run('EQL0D_initialize.m'); %%% Initialize or restart from .mat file
    
    %% Outer Loop / Cycles
    while(~SYS.stopOuter)
        SYS.ouCntr=SYS.ouCntr+1; SYS.inCntr=0;
        if(strcmp(OPT.iterMode,'equilibrium'))
            SYS.tStep(end+1)=OPT.cycleLength(SYS.ouCntr)*24.0*3600.0;
        elseif(strcmp(OPT.iterMode,'steps'))
            SYS.tStep(end+1)=OPT.cycleLength(SYS.ouCntr)*24.0*3600.0/OPT.nSteps(SYS.ouCntr);
        end
        SYS.stopInner=false;  SYS.PCC.corrector=false;
        SYS.oldFIMA=[MAT.FIMA]; SYS.oldN = [MAT.N]; %%% Save previous cycle data
        if(~SYS.debugMode)
            %%% Adapt depletion time/burnup in Serpent
            modifyInput(SYS.Casename,'dep',OPT.cycleLength(SYS.ouCntr));
            runSerpent(OPT,SYS);
        end
        if(~SYS.debugMode||SYS.ouCntr==1)
            [MAT,SYS] = loadSerpentData(MAT,SYS); %%% Read Serpent outputs
            MAT=updateRates(MAT,SYS);
            if(~SYS.debugMode)
                saveFiles(MAT,OPT,SYS);  %%% Move files to folder
                printK(SYS,'AB','C','Serpent');
            end
        end
        SYS = computeK(MAT,SYS); %%% Compute k-eff
        printK(SYS,'AB','C','EQL0D');
        if(~isempty(SYS.IDX.contStr))
            SYS = createRepMatrices(MAT,REP,SYS);
        end
        SYS = buildSystemMatrix(SYS);
        save([SYS.Casename '.mat']);
        while(~SYS.stopInner) %%% Inner loop
            SYS.inCntr=SYS.inCntr+1; SYS.prevFIMA=[MAT(SYS.IDX.contMat).FIMA];
            SYS.prevN.BOC=[MAT.N]; %%% Store compositions from previous loop
            
            [MAT,SYS]=burnCycle(MAT,OPT,REP,SYS);
            
            if(OPT.PCC&&~SYS.debugMode) %Corrector step
                SYS.PCC.corrector=true; SYS.PCC.nSteps=OPT.nSteps(SYS.ouCntr);
                runSerpent(OPT,SYS); %%% Run Serpent/Read outputs
                [MAT,SYS] = loadSerpentData(MAT,SYS);
                MAT=updateRates(MAT,SYS);
                if(~SYS.debugMode)
                    saveFiles(MAT,OPT,SYS);
                    printK(SYS,'AB','P','Serpent');
                end
                SYS = computeK(MAT,SYS);
                printK(SYS,'AB','P','EQL0D');
                
                if(~isempty(SYS.IDX.contStr))
                    SYS = createRepMatrices(MAT,REP,SYS);
                end
                SYS = buildSystemMatrix(SYS);
                
                while(SYS.inCntr<=OPT.nSteps(SYS.ouCntr)) %%% Burn system
                    [MAT,SYS]=burnCycle(MAT,OPT,REP,SYS);
                    SYS.inCntr=SYS.inCntr+1;
                end
                SYS.PCC.corrector=false; SYS.stopInner=true;
            else
                SYS=testConvergence(MAT,OPT,SYS,'inner');
            end
            %save([SYS.Casename '.mat']);
        end
        if(OPT.printCycles&&~OPT.printSteps)
            for i=SYS.IDX.contMat
                MAT(i).printMaterial(SYS,'EoC'); %%% Print compositions at EoC
            end
        end
        SYS=testConvergence(MAT,OPT,SYS,'outer'); %%% stop outer loop?
    end
    
    %% Last step / Equilibrium results
    if(~SYS.debugMode)
        runSerpent(OPT,SYS);
    end
    [MAT,SYS] = loadSerpentData(MAT,SYS);
    switch OPT.iterMode
        case 'equilibrium'
            for i=SYS.IDX.contMat
                MAT(i).printMaterial(SYS,'EQL_AB');
                mat=MAT(i); mat.N=SYS.prevN.EOC(:,i);
                mat.printMaterial(SYS,'EQL_BB');
            end
            neutronBalance(MAT,SYS)
        case 'steps'
            if(OPT.printSteps)
                for i=SYS.IDX.contMat
                    MAT(i).printMaterial(SYS,'EoC');
                end
            end
    		if(~SYS.debugMode)
        		saveFiles(MAT,OPT,SYS);
    		end
    end
    save([SYS.Casename '.mat']);
    fprintf(SYS.FID.log,'%s\n','**** EQL0D **** Procedure finished.');
    for file=fields(SYS.FID)
        fclose(SYS.FID.(file{1}));
    end
    if(OPT.writeMail)
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
