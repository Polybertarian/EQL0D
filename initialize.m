function [MAT,OPT,REP,SYS] = initialize(SYS)
%[MAT,OPT,REP,SYS] = INITIALIZE(SYS) regroups all steps to parse the user input and initialize
%the initial variables of EQL0D

stdFields={'nCores','debugMode','verboseMode','printAndQuit','restartCalc','restartCalc','resetCounters'};
for i=find(~isfield(SYS,stdFields))
    SYS.(stdFields{i})=false;
end

if(SYS.restartCalc)
    load([SYS.Casename '.mat']);
    [SYS.FID.log,errmsg]=fopen([SYS.Casename '.log'],'at');
    if(~isempty(errmsg))
        error(errmsg)
    end
    [SYS.FID.keff,errmsg]=fopen('keff.txt','at');

    if(OPT.reactControl)
        switch OPT.REA.mode
            case {'replace','addMass'}
                SYS.FID.react=fopen('reactivity.txt','at');
            case 'addVolume'
                [SYS.FID.volume,errmsg]=fopen('volume.txt','at');
        end
    end
    if(~isempty(errmsg))
        error(errmsg)
    end
else
    run('defaultConfig.m'); %%% Default config
    if(exist([SYS.Casename '.m'],'file')==2)
        run([SYS.Casename '.m']);
    else
        error(['Cannot find configuration file ' SYS.Casename '.m !'])
    end
    [SYS.FID.log,errmsg]=fopen([SYS.Casename '.log'],'wt');
    if(~isempty(errmsg))
        error(errmsg)
    end
    if(SYS.resetCounters)
        tmp=load([SYS.Casename '.mat'],'MAT');
        for i=find([MAT.isCont]&~[MAT.isStr])
            MAT(i).N(:,end+1)=tmp.MAT(i).N(:,end);
        end
    end

    %%% load isotope list and properties
    if(exist('DAT','var')~=1)
        warning('Library undefined! Loading default...')
        load([OPT.defaultDataLibrary '.mat'],DAT)
    end
    if(~exist('MAT','var'))
        error('Error: MAT vector undefined!')
    elseif(isempty(MAT))
        error('Error: MAT vector empty!')
    end
    if(~exist('REP','var'))
        warning('No reprocessing streams defined!')
    end

    %%% Initialization of various variables
    SYS.ouCntr=0;
    SYS.inCntr=0;
    SYS.stopOuter=false;
    SYS.nowTime=0.0;
    SYS.tStep=[];
    [SYS.RR.notInMat{1}]=deal(struct('fiss',[],'capt',[],'n2n',[],'n3n',[]));
    [SYS.RR.inMat{1}]=deal(struct('fiss',zeros(size(DAT.ZAI0)),'capt',zeros(size(DAT.ZAI0)),...
        'n2n',zeros(size(DAT.ZAI0)),'n3n',zeros(size(DAT.ZAI0))));
    SYS.keff=[];
    SYS.kinf=[];
    SYS.RR.NU=cell(2,1);
    SYS.RR.LEAK=cell(2,1);
    SYS.MTX.defaultDecay=sparse(DAT.decayMatrix(isProduced(DAT.ZAI0),isProduced(DAT.ZAI0)));

    %%% Write material compositions for Serpent
    if(SYS.printAndQuit)
        for i=find([MAT.isBurned])
            MAT(i).write('');
            MAT(i).printMaterial(SYS,'EoC');
        end
        save([SYS.Casename '.mat']);
        exit(0)
    else
        %%% Parse all data
        [MAT,OPT,REP,SYS] = parseInput(MAT,OPT,REP,SYS);
        for i=SYS.IDX.burnMat
            MAT(i).write(OPT.matWriteStyle);
        end
        for i=[SYS.IDX.burnMat SYS.IDX.strMat]
            if(OPT.printSteps)
                MAT(i).printMaterial(SYS,'AB');
            elseif(OPT.printCycles&&~OPT.printSteps)
                MAT(i).printMaterial(SYS,'EoC');
            end
        end
    end

    %%% Neutron balance outputs
    [SYS.FID.keff,~]=fopen('keff.txt','w');
    fprintf(SYS.FID.keff,'%-7s %-3s %-5s %-4s %-3s %-12s %-9s %-9s\n',...
        'Source','PCC','Cycle','Step','Rep','Time','k-inf','k-eff');

    %%% Reactivity control outputs
    if(OPT.reactControl)
        switch OPT.REA.mode
            case {'addMass','replace'}
                SYS.FID.react=fopen('reactivity.txt','w');
                name{1}=MAT(SYS.IDX.targetMat).name;
                if(~isempty(SYS.IDX.feedMat))
                    name{2}=MAT(SYS.IDX.feedMat).name;
                else
                    name{2}=[];
                end
                fprintf(SYS.FID.react,'%-22s','Time');
                fprintf(SYS.FID.react,['%-' num2str(12*numel(SYS.IDX.targetNucUp)) 's'],name{1});
                if(strcmp(OPT.REA.mode,'addMass')&&~isempty(SYS.IDX.feedMat))
                    fprintf(SYS.FID.react,['%-' num2str(13*numel(SYS.IDX.targetNucUp)) 's'],name{2});
                elseif(strcmp(OPT.REA.mode,'replace'))
                    fprintf(SYS.FID.react,['%-' num2str(13*numel(SYS.IDX.targetNucRepl)) 's'],name{1});
                    if(~isempty(SYS.IDX.feedMat))
                        fprintf(SYS.FID.react,['%-' num2str(13*numel(SYS.IDX.targetNucUp)) 's'],name{2});
                        fprintf(SYS.FID.react,['%-' num2str(13*numel(SYS.IDX.targetNucRepl)) 's'],name{2});
                    end
                end
                fprintf(SYS.FID.react,'\n');
                fprintf(SYS.FID.react,'%-7s%-6s%-9s','Cycle','Step','EFPD');
                names{1}=MAT(SYS.IDX.targetMat).nuclideName(SYS.IDX.targetNucUp);
                fprintf(SYS.FID.react,repmat('%-13s',1,numel(names{1})),names{1}{:});
                if(strcmp(OPT.REA.mode,'addMass')&&~isempty(SYS.IDX.feedMat))
                    names{2}=MAT(SYS.IDX.feedMat).nuclideName(SYS.IDX.feedNucUp);
                    fprintf(SYS.FID.react,repmat('%-13s',1,numel(names{2})),names{2}{:});
                elseif(strcmp(OPT.REA.mode,'replace'))
                    names{2}=MAT(SYS.IDX.targetMat).nuclideName(SYS.IDX.targetNucDo);
                    fprintf(SYS.FID.react,repmat('%-13s',1,numel(names{2})),names{2}{:});
                    if(~isempty(SYS.IDX.feedMat))
                        names{3}=MAT(SYS.IDX.feedMat).nuclideName(SYS.IDX.feedNucUp);
                        fprintf(SYS.FID.react,repmat('%-13s',1,numel(names{3})),names{3}{:});
                        names{4}=MAT(SYS.IDX.feedMat).nuclideName(SYS.IDX.feedNucDo);
                        fprintf(SYS.FID.react,repmat('%-13s',1,numel(names{4})),names{4}{:});
                    end
                end
                fprintf(SYS.FID.react,'\n');
            case 'addVolume'
                [SYS.FID.volume,~]=fopen('volume.txt','w');
                fprintf(SYS.FID.volume,'%-7s%-6s%-9s%-12s\n',...
                    'Cycle','Step','EFPD','Add. Vol. [cm^3]');
        end
    end

    %%% Check presence of necessary input parameters in Serpent file
    [~,isAbsent]=unix(['grep -c "set depmtx 1" ' SYS.Casename]);
    if(str2double(isAbsent)==0)
        [fID,~]=fopen(SYS.Casename,'a');
        fprintf(fID,'\n%s\n','set depmtx 1 1');
        fclose(fID);
    end
    [~,isAbsent]=unix(['grep -c "set arr 1" ' SYS.Casename]);
    if(str2double(isAbsent)==0)
        [fID,~]=fopen(SYS.Casename,'a');
        fprintf(fID,'\n%s\n','set arr 1 0');
        fclose(fID);
    end
    [~,isAbsent]=unix(['grep -c "det intFlux" ' SYS.Casename]);
    if(str2double(isAbsent)==0)
        [fID,~]=fopen(SYS.Casename,'a');
        fprintf(fID,'\n%s','det intFlux ');
        for i=SYS.IDX.fluxMat
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
