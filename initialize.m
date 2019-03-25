function [MAT,OPT,REP,SYS] = initialize(SYS)
    %[MAT,OPT,REP,SYS] = INITIALIZE(SYS) regroups all steps to parse the user input and initialize
    %the initial variables of EQL0D

    stdFields={'nCores','debugMode','verboseMode','printAndQuit','restartCalc','restartCalc','resetCounters'};
    for i=find(~isfield(SYS,stdFields))
        SYS.(stdFields{i})=false;
    end

    if SYS.restartCalc
        load([SYS.Casename '.mat']);
        [SYS.FID.log,errmsg]=fopen([SYS.Casename '.log'],'at');
        if ~isempty(errmsg)
            error(errmsg)
        end
        [SYS.FID.keff,errmsg]=fopen('keff.txt','at');

        if OPT.reactControl
            switch OPT.REA.mode
            case {'replace','addMass'}
                SYS.FID.react=fopen('reactivity.txt','at');
            case 'addVolume'
                [SYS.FID.volume,errmsg]=fopen('volume.txt','at');
            end
        end
        if ~isempty(errmsg)
            error(errmsg)
        end
    else
        run('defaultConfig.m'); % Default config
        run('userPrefs.m');
        if exist([SYS.Casename '.m'],'file')==2
            run([SYS.Casename '.m']);
        else
            error(['Cannot find configuration file ' SYS.Casename '.m !'])
        end
        [SYS.FID.log,errmsg]=fopen([SYS.Casename '.log'],'wt');
        if ~isempty(errmsg)
            error(errmsg)
        end
        if SYS.resetCounters
            tmp=load([SYS.Casename '.mat'],'MAT');
            for i=find([MAT.isCont]&~[MAT.isStr])
                MAT(i).N(:,end+1)=tmp.MAT(i).N(:,end);
            end
        end

        if exist('DAT','var')~=1 % load isotope list and properties
            warning('Library undefined! Loading default...')
            load([OPT.defaultDataLibrary '.mat'],DAT)
        end
        SYS.nuclearDataLibrary=DAT.libraryName;
        if ~exist('MAT','var')
            error('Error: MAT vector undefined!')
        elseif(isempty(MAT))
            error('Error: MAT vector empty!')
        end
        if ~exist('REP','var')
            warning('No reprocessing streams defined!')
        end

        SYS.ouCntr=0; SYS.inCntr=0; % Initialization of various variables
        SYS.stopOuter=false; SYS.nowTime=0.0; SYS.tStep=[];
        SYS.KEFF.EQL0D=[];  SYS.KEFF.Serpent=[]; SYS.KINF.EQL0D=[];  SYS.KINF.Serpent=[];
        SYS.MTX.defaultDecay=sparse(DAT.decayMatrix(isProduced(DAT.libraryName,DAT.ZAI0),isProduced(DAT.libraryName,DAT.ZAI0)));

        if SYS.printAndQuit         % Write material compositions for Serpent
            for i=find([MAT.isBurned])
                MAT(i).write('');
                MAT(i).printMaterial(SYS,'EoC');
            end
            save([SYS.Casename '.mat']);
            exit(0)
        else
            %%% Parse all data
            [MAT,OPT,REP,SYS] = parseInput(MAT,OPT,REP,SYS);
            for i=SYS.IDX.MAT.burn
                MAT(i).write(OPT.matWriteStyle);
            end
            for i=[SYS.IDX.MAT.burn SYS.IDX.MAT.decay]
                if OPT.printSteps
                    MAT(i).printMaterial(SYS,'AB');
                elseif OPT.printCycles&&~OPT.printSteps
                    MAT(i).printMaterial(SYS,'EoC');
                end
            end
        end
        SYS.RR.NU=zeros(3,numel(SYS.IDX.MAT.inFlux));
        SYS.RR.LEAK=zeros(3,1);
        SYS.RR.devFiss=zeros(3,numel(SYS.IDX.MAT.inFlux));
        SYS.RR.devCapt=zeros(3,numel(SYS.IDX.MAT.inFlux));
        [SYS.RR.notInMat{1}]=deal(struct('fiss',[],'capt',[],'n2n',[],'n3n',[]));
        SYS.RR.notInMat{2}=SYS.RR.notInMat{1}; SYS.RR.notInMat{3}=SYS.RR.notInMat{2};
        [SYS.RR.inMat{1}]=deal(struct('fiss',zeros(size(MAT(1).ZAI)),'capt',zeros(size(MAT(1).ZAI)),...
        'n2n',zeros(size(MAT(1).ZAI)),'n3n',zeros(size(MAT(1).ZAI))));
        SYS.RR.inMat{2}=SYS.RR.inMat{1}; SYS.RR.inMat{3}=SYS.RR.inMat{1};

        %%% Neutron balance outputs
        [SYS.FID.keff,~]=fopen('keff.txt','w');
        fprintf(SYS.FID.keff,'%-7s %-3s %-5s %-4s %-3s %-12s %-9s %-9s\n','Source','PCC','Cycle','Step','Rep','Time','k-inf','k-eff');

        %%% Reactivity control outputs
        if OPT.reactControl
            switch OPT.REA.mode
            case {'addMass','replace'}
                SYS.FID.react=fopen('reactivity.txt','w');
                name{1}=MAT(SYS.IDX.REA.target).name;
                if ~isempty(SYS.IDX.REA.feed)
                    name{2}=MAT(SYS.IDX.REA.feed).name;
                else
                    name{2}=[];
                end
                fprintf(SYS.FID.react,'%-22s','Time');
                if strcmp(OPT.REA.mode,'addMass')&&~isempty(SYS.IDX.REA.feed)
                    elementsTarget=unique([SYS.IDX.REA.targetNucUp;SYS.IDX.REA.targetNucDo]);
                    fprintf(SYS.FID.react,['%-' num2str(13*numel(elementsTarget)) 's'],name{1});
                elseif strcmp(OPT.REA.mode,'replace')
                    elementsTarget=unique([SYS.IDX.REA.targetNucRepl;SYS.IDX.REA.targetNucUp;SYS.IDX.REA.targetNucDo]);
                    fprintf(SYS.FID.react,['%-' num2str(13*numel(elementsTarget)) 's'],name{1});
                end
                if ~isempty(SYS.IDX.REA.feed)
                    fprintf(SYS.FID.react,['%-' num2str(13*numel(elementsTarget)) 's'],name{2});
                end
                fprintf(SYS.FID.react,'\n');
                fprintf(SYS.FID.react,'%-7s%-6s%-9s','Cycle','Step','EFPD');
                names{1}=MAT(SYS.IDX.REA.target).nuclideName(elementsTarget);
                fprintf(SYS.FID.react,repmat('%-13s',1,numel(names{1})),names{1}{:});
                if ~isempty(SYS.IDX.REA.feed)
                    fprintf(SYS.FID.react,repmat('%-13s',1,numel(names{1})),names{1}{:});
                end
                fprintf(SYS.FID.react,'\n');
            case 'addVolume'
                [SYS.FID.volume,~]=fopen('volume.txt','w');
                fprintf(SYS.FID.volume,'%-7s%-6s%-9s%-12s\n','Cycle','Step','EFPD','Add. Vol. [cm^3]');
            end
        end

        %%% Check presence of necessary input parameters in Serpent file
        %[~,isAbsent]=unix(['grep -c "set depmtx 1" ' SYS.Casename]);
        %if(str2double(isAbsent)==0)
        %  [fID,~]=fopen(SYS.Casename,'a');
        %  fprintf(fID,'\n%s\n','set depmtx 1 1');
        %  fclose(fID);
        %end
        [~,isAbsent]=unix(['grep -c "set arr 1" ' SYS.Casename]);
        if str2double(isAbsent)==0
            [fID,~]=fopen(SYS.Casename,'a');
            fprintf(fID,'\n%s\n','set arr 1 0');
            fclose(fID);
        end
        [~,isAbsent]=unix(['grep -c "det intFlux" ' SYS.Casename]);
        if str2double(isAbsent)==0
            [fID,~]=fopen(SYS.Casename,'a');
            fprintf(fID,'\n%s','det intFlux ');
            for i=SYS.IDX.MAT.inFlux
                fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
            end
            fclose(fID);
        end
    end
    clearvars -except OPT SYS DAT MAT REP
    save([SYS.Casename '.mat']);
    fprintf(SYS.FID.log,'%s\n','**** EQL0D **** Procedure initialized.');

    return
end
