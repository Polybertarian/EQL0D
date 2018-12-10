function [MAT,OPT,REP,SYS] = initialize(SYS)
%[MAT,OPT,REP,SYS] = INITIALIZE(SYS) regroups all steps to parse the user input and initialize
%the initial variables of EQL0D
global FID

try
    if SYS.restartCalc
        load([SYS.Casename '.mat']);
        [FID.log,errmsg]=fopen([Casename '.log'],'at');
        openLogs(true,SYS.REA); %re-open logs
        fprintf(FID.log,'%s\n','**** EQL0D **** Procedure re-started.');
    else
        stdFields={'nCores','debugMode','verboseMode','printAndQuit','restartCalc','restartCalc','resetCounters'};
        for i=find(~isfield(SYS,stdFields))
            SYS.(stdFields{i})=false;
        end
        [FID.log,errmsg]=fopen([Casename '.log'],'wt');
        fprintf(FID.log,'%s\n','**** EQL0D **** Initializing procedure...');
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
        
        if exist('DAT','var')~=1 % load isotope list and properties
            fprintf(FID.log,'%s\n',['**** EQL0D **** Library undefined! Loading default ''' OPT.defaultDataLibrary '''']);
            load([OPT.defaultDataLibrary '.mat'],DAT)
        end
        SYS.nuclearDataLibrary=DAT.libraryName;
        if ~exist('MAT','var')
            error('EQLOD:MATVectorUndefined','Error: MAT vector undefined!')
        elseif isempty(MAT)
            error('EQLOD:MATVectorEmpty','Error: MAT vector empty!')
        end
        if ~exist('REP','var')
            fprintf(FID.log,'%s\n','**** EQL0D **** Warning: No reprocessing streams defined!');
            REP=[];
        elseif isempty(REP)
            fprintf(FID.log,'%s\n','**** EQL0D **** Warning: REP vector empty. No reprocessing streams defined!');
        end
        
        SYS.ouCntr=0; SYS.inCntr=0; SYS.stopOuter=false;  SYS.nowTime=0.0; SYS.tStep=[];
        tmpVec=zeros(size(MAT(1).ZAI));
        SYS.RR=struct('inMat',deal(struct('fiss',tmpVec,'capt',tmpVec,'n2n',tmpVec,'n3n',tmpVec)),...
            'notInMat',deal(struct('fiss',[],'capt',[],'n2n',[],'n3n',[])),...
            'NU',zeros(3,1),'LEAK',zeros(3,1),'devFiss',zeros(3,1),'devCapt',zeros(3,1),...
            'KEFF',struct('EQL0D',[],'Serpent',[]),'KINF',struct('EQL0D',[],'Serpent',[]));
        SYS.RR.notInMat{2}=SYS.RR.notInMat{1};   SYS.RR.notInMat{3}=SYS.RR.notInMat{1};
        SYS.RR.inMat{2}=SYS.RR.inMat{1};  SYS.RR.inMat{3}=SYS.RR.inMat{1};
        
        if SYS.resetCounters
            tmp=load([SYS.Casename '.mat'],'MAT');
            for i=find([MAT.isCont]&~[MAT.isStr])
                MAT(i).N(:,end+1)=tmp.MAT(i).N(:,end);
            end
        end
        
        [MAT,OPT,REP,SYS] = parseInput(MAT,OPT,REP,SYS);  % Parse all data
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
    fprintf(FID.log,'%s\n','**** EQL0D **** Procedure initialized.');
catch
    fprintf(FID.log,'%s\n','**** EQL0D **** Error during initialization!');
    exit(1);
end

return
end
