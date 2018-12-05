function [] = openLogs(Casename,restartCalc,REA)
%OPENLOGS Opens log files for EQL0D
global FID

if(restartCalc)
  [FID.log,errmsg]=fopen([Casename '.log'],'at');
  [FID.keff,errmsg]=fopen('keff.txt','at');
  if(REA.reactControl)
    switch REA.mode
      case {'replace','addMass'}
        FID.react=fopen('reactivity.txt','at');
      case 'addVolume'
        [FID.volume,errmsg]=fopen('volume.txt','at');
    end
  end
else
  [FID.log,errmsg]=fopen([Casename '.log'],'wt');
  [FID.keff,errmsg]=fopen('keff.txt','w');
  fprintf(FID.keff,'%-7s %-3s %-5s %-4s %-3s %-12s %-9s %-9s\n',...
    'Source','PCC','Cycle','Step','Rep','Time','k-inf','k-eff');
  if(REA.reactControl)
    switch REA.mode
      case {'addMass','replace'}
        FID.react=fopen('reactivity.txt','w');
        name{1}=REA.targetMat;
        if(~isempty(REA.feedMat))
          name{2}=REA.feedMat;
        else
          name{2}=[];
        end
        fprintf(FID.react,'%-22s','Time');
        elementsTarget=unique([REA.replNuclides;REA.upNuclides;REA.downNuclides]);
        fprintf(FID.react,['%-' num2str(13*numel(elementsTarget)) 's'],name{1});
        if(~isempty(REA.feedMat))
          fprintf(FID.react,['%-' num2str(13*numel(elementsTarget)) 's'],name{2});
        end
        fprintf(FID.react,'\n');
        fprintf(FID.react,'%-7s%-6s%-9s','Cycle','Step','EFPD');
        names{1}=ZAI2Name(elementsTarget);
        fprintf(FID.react,repmat('%-13s',1,numel(names{1})),names{1}{:});
        if(~isempty(REA.feedMat))
          fprintf(FID.react,repmat('%-13s',1,numel(names{1})),names{1}{:});
        end
        fprintf(FID.react,'\n');
      case 'addVolume'
        [FID.volume,~]=fopen('volume.txt','w');
        fprintf(FID.volume,'%-7s%-6s%-9s%-12s\n','Cycle','Step','EFPD','Add. Vol. [cm^3]');
    end
  end
end
end

