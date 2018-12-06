function [MAT,OPT,REP,SYS] = initialize(SYS)
%[MAT,OPT,REP,SYS] = INITIALIZE(SYS) regroups all steps to parse the user input and initialize
%the initial variables of EQL0D

stdFields={'nCores','debugMode','verboseMode','printAndQuit','restartCalc','restartCalc','resetCounters'};
for i=find(~isfield(SYS,stdFields))
    SYS.(stdFields{i})=false;
end

if(SYS.restartCalc)
    load([SYS.Casename '.mat']);
    openLogs(SYS.Casename,true,SYS.REA);
else
    run('defaultConfig.m'); % Default config
    if(exist([SYS.Casename '.m'],'file')==2)
        run([SYS.Casename '.m']);
    else
        error(['Cannot find configuration file ' SYS.Casename '.m !'])
    end
    
    if(SYS.resetCounters)
        tmp=load([SYS.Casename '.mat'],'MAT');
        for i=find([MAT.isCont]&~[MAT.isStr])
            MAT(i).N(:,end+1)=tmp.MAT(i).N(:,end);
        end
    end
    
    if(exist('DAT','var')~=1) % load isotope list and properties
        warning('Library undefined! Loading default...')
        load([OPT.defaultDataLibrary '.mat'],DAT)
    end
    SYS.nuclearDataLibrary=DAT.libraryName;
    if(~exist('MAT','var'))
        error('Error: MAT vector undefined!')
    elseif(isempty(MAT))
        error('Error: MAT vector empty!')
    end
    if(~exist('REP','var'))
        warning('No reprocessing streams defined!')
    end
    
    SYS.ouCntr=0; SYS.inCntr=0; SYS.stopOuter=false; % Initialization of various variables
    SYS.nowTime=0.0; SYS.tStep=[];
    [SYS.RR.notInMat{1}]=deal(struct('fiss',[],'capt',[],'n2n',[],'n3n',[]));
    SYS.RR.notInMat{2}=SYS.RR.notInMat{1};   SYS.RR.notInMat{3}=SYS.RR.notInMat{1};
    tmpVec=zeros(size(MAT(1).ZAI));
    [SYS.RR.inMat{1}]=deal(struct('fiss',tmpVec,'capt',tmpVec,'n2n',tmpVec,'n3n',tmpVec));
    SYS.RR.inMat{2}=SYS.RR.inMat{1};SYS.RR.inMat{3}=SYS.RR.inMat{1};
    SYS.RR.NU=zeros(3,1);SYS.RR.LEAK=zeros(3,1);
    SYS.RR.devFiss=zeros(3,1); SYS.RR.devCapt=zeros(3,1);
    SYS.KEFF.EQL0D=[];  SYS.KEFF.Serpent=[];
    SYS.KINF.EQL0D=[];  SYS.KINF.Serpent=[];
    if OPT.PCC
        OPT.renormalize=false;
    end
   
    if(SYS.printAndQuit) % Write material compositions for Serpent
        for i=find([MAT.isBurned])
            MAT(i).write('');
            MAT(i).printMaterial(SYS,'EoC');
        end
        save([SYS.Casename '.mat']);
        exit(0)
    else
        if(exist('REP','var')==0)
            REP=[];
        end
        [MAT,OPT,REP,SYS] = parseInput(MAT,OPT,REP,SYS);  % Parse all data
        for i=SYS.IDX.MAT.burn
            MAT(i).write(OPT.matWriteStyle);
        end
        for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
            if(OPT.printSteps)
                MAT(i).printMaterial(SYS,'AB');
            elseif(OPT.printCycles&&~OPT.printSteps)
                MAT(i).printMaterial(SYS,'EoC');
            end
        end
    end
    
    openLogs(SYS.Casename,false,SYS.REA);
    
    isAbsent=[]; % Check presence of necessary input parameters in Serpent file
    %[~,isAbsent{end+1}]=unix(['grep -c "set depmtx 1" ' SYS.Casename]);
    [~,isAbsent{end+1}]=unix(['grep -c "set arr 1" ' SYS.Casename]);
    [~,isAbsent{end+1}]=unix(['grep -c "det intFlux" ' SYS.Casename]);
    [~,isAbsent{end+1}]=unix(['grep -c "det intFiss" ' SYS.Casename]);
    [~,isAbsent{end+1}]=unix(['grep -c "det intProd" ' SYS.Casename]);
    [~,isAbsent{end+1}]=unix(['grep -c "det intCapt" ' SYS.Casename]);
    isAbsent=logical(str2double(isAbsent));
    
    if(any(isAbsent))
        [fID,~]=fopen(SYS.Casename,'a');
        % if(isAbsent(1))
        % fprintf(fID,'\n%s\n','set depmtx 1 1'); 
        % end
        if(isAbsent(1))
            fprintf(fID,'\n%s\n','set arr 1 0');
        end
        if(isAbsent(2))
            fprintf(fID,'\n%s','det intFlux ');
            for i=SYS.IDX.MAT.inFlux
                fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
            end
        end
        if(isAbsent(3))
            fprintf(fID,'\n%s','det intFiss dr -6 void ');
            for i=SYS.IDX.MAT.inFlux
                fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
            end
        end
        if(isAbsent(4))
            fprintf(fID,'\n%s','det intProd dr -7 void ');
            for i=SYS.IDX.MAT.inFlux
                fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
            end
        end
        if(isAbsent(5))
            fprintf(fID,'\n%s','det intCapt dr -2 void ');
            for i=SYS.IDX.MAT.inFlux
                fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
            end
        end
        fclose(fID);
    end
end
clearvars -except OPT SYS DAT MAT REP
save([SYS.Casename '.mat']);
fprintf(SYS.FID.log,'%s\n','**** EQL0D **** Procedure initialized.');

return
end
