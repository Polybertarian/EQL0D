function [MAT,OPT,REP,SYS] = initialize(SYS)
%[MAT,OPT,REP,SYS] = INITIALIZE(SYS) regroups all steps to parse the user input and initialize
%the initial variables of EQL0D

global FID

stdFields={'nCores','debugMode','verboseMode','printAndQuit','restartCalc','restartCalc','resetCounters'};
for i=find(~isfield(SYS,stdFields))
  SYS.(stdFields{i})=false;
end

if(SYS.restartCalc)
  load([SYS.Casename '.mat']);
  openLogs(SYS.Casename,true,SYS.REA)
else
  run('defaultConfig.m'); %%% Default config
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
  
  %%% load isotope list and properties
  if(exist('DAT','var')~=1)
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
  
  %%% Initialization of various variables
  SYS.ouCntr=0;  SYS.inCntr=0;
  SYS.stopOuter=false;  SYS.nowTime=0.0;  SYS.tStep=[];
  [SYS.RR.notInMat{1}]=deal(struct('fiss',[],'capt',[],'n2n',[],'n3n',[]));
  SYS.RR.notInMat{2}=SYS.RR.notInMat{1};
  tmpVec=zeros(size(MAT(1).ZAI));
  [SYS.RR.inMat{1}]=deal(struct('fiss',tmpVec,'capt',tmpVec,'n2n',tmpVec,'n3n',tmpVec));
  SYS.RR.inMat{2}=SYS.RR.inMat{1};
  SYS.KEFF.EQL0D=[];  SYS.KEFF.Serpent=[];
  SYS.KINF.EQL0D=[];  SYS.KINF.Serpent=[];
  SYS.RR.NU=cell(2,1);
  SYS.RR.LEAK=cell(2,1);
  if OPT.PCC
    OPT.renormalize=false;
  end
  
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
    openLogs(SYS.Casename,false,SYS.REA)
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
  
  %%% Check presence of necessary input parameters in Serpent file
  %[~,isAbsent]=unix(['grep -c "set depmtx 1" ' SYS.Casename]);
  %if(str2double(isAbsent)==0)
  %  [fID,~]=fopen(SYS.Casename,'a');
  %  fprintf(fID,'\n%s\n','set depmtx 1 1');
  %  fclose(fID);
  %end
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
    for i=SYS.IDX.MAT.inFlux
      fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
    end
    fclose(fID);
  end
end
clearvars -except OPT SYS DAT MAT REP FID
save([SYS.Casename '.mat']);
fprintf(FID.log,'%s\n','**** EQL0D **** Procedure initialized.');

return
end
