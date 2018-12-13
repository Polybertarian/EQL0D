function [MAT,OPT,REP,SYS] = initialize(SYS)
%[MAT,OPT,REP,SYS] = INITIALIZE(SYS) regroups all steps to parse the user input and initialize
%the initial variables of EQL0D
global FID

try
    stdFields={'nCores','debugMode','verboseMode','printAndQuit','restartCalc','restartCalc','resetCounters'};
    for i=find(~isfield(SYS,stdFields))
        SYS.(stdFields{i})=false;
    end
    if SYS.restartCalc
        load([SYS.Casename '.mat']);
        [FID.log,errmsg]=fopen([SYS.Casename '.log'],'at');
        openLogs(true,SYS.REA); %re-open logs
        fprintf('%s\n','**** EQL0D **** Procedure re-started.');
    else
        [FID.log,errmsg]=fopen([SYS.Casename '.log'],'wt');
        fprintf('%s\n','**** EQL0D **** Initializing procedure...');
        run('defaultConfig.m'); % Default config
        run('userPrefs.m'); % User preferences
        if exist([SYS.Casename '.m'],'file')==2
            run([SYS.Casename '.m']);
        else
            error('EQL0D:InputNotFound',['Cannot find EQL0D configuration file ' SYS.Casename '.m !'])
        end
        
        if SYS.printAndQuit
            for i=find([MAT.isBurned])
                MAT(i).write(''); % Write material compositions for Serpent
                MAT(i).printMaterial(0,0,0,'EoC')
            end
            save([SYS.Casename '.mat']);
            exit(0)
        end
        
        if exist('DAT','var')~=1 % Load isotope list and properties
            fprintf('%s\n',['**** EQL0D **** Library undefined! Loading default ''' OPT.defaultDataLibrary '''']);
            load([OPT.defaultDataLibrary '.mat'],DAT)
        end
        SYS.nuclearDataLibrary=DAT.libraryName;
        if ~exist('MAT','var')
            error('EQLOD:MATVectorUndefined','Error: MAT vector undefined!')
        elseif isempty(MAT)
            error('EQLOD:MATVectorEmpty','Error: MAT vector empty!')
        end
        if ~exist('REP','var')
            fprintf('%s\n','**** EQL0D **** Warning: No reprocessing streams defined!');
            REP=[];
        elseif isempty(REP)
            fprintf('%s\n','**** EQL0D **** Warning: REP vector empty. No reprocessing streams defined!');
        end
        
        SYS.ouCntr=0; SYS.inCntr=0; SYS.stopOuter=false;  SYS.nowTime=0.0; SYS.tStep=[];
        
        if SYS.resetCounters
            tmp=load([SYS.Casename '.mat'],'MAT');
            for i=find([MAT.isCont]&~[MAT.isStr])
                MAT(i).N(:,end+1)=tmp.MAT(i).N(:,end);
            end
        end
        
        [MAT,OPT,REP,SYS] = parseInput(MAT,OPT,REP,SYS);  % Parse all data
        tmpVec=zeros(size(MAT(1).ZAI)); tmpZero=zeros(1,numel(SYS.IDX.MAT.inFlux));
        SYS.RR=struct('inMat',struct('fiss',tmpVec,'capt',tmpVec,'n2n',tmpVec,'n3n',tmpVec),...
            'notInMat',struct('fiss',zeros(1),'capt',zeros(1),'n2n',zeros(1),'n3n',zeros(1)),...
            'NU',tmpZero,'LEAK',zeros(1),'devFiss',tmpZero,'devCapt',tmpZero);
        SYS.RR(2)=SYS.RR(1);SYS.RR(3)=SYS.RR(2);
        SYS.KEFF=struct('EQL0D',[],'Serpent',[]);SYS.KINF=struct('EQL0D',[],'Serpent',[]);
        
        openLogs(false,SYS.REA);
        for i=SYS.IDX.MAT.burn
            MAT(i).write(OPT.matWriteStyle);
        end
        for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
            if OPT.printSteps
                MAT(i).printMaterial(SYS.ouCntr,SYS.inCntr,SYS.nowTime,'AB');
            elseif OPT.printCycles&&~OPT.printSteps
                MAT(i).printMaterial(SYS.ouCntr,SYS.inCntr,SYS.nowTime,'EoC');
            end
        end
        modifySerpentInput(SYS.Casename,'det',[])
    end
    clearvars -except OPT SYS DAT MAT REP FID
    save([SYS.Casename '.mat']);
    fprintf('%s\n','**** EQL0D **** Procedure initialized.');
catch exception
    fprintf('%s\n','**** EQL0D **** Error during initialization!');
    rethrow(exception)
end

return
end
