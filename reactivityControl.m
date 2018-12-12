function [MAT,SYS] = reactivityControl(MAT,SYS)
%REACTIVITYCONTROL(MAT,OPT,SYS) Controls reactivity by changing the concentrations of
%selected nuclides in selected materials of the MAT vector according to
%user options given in the OPT.REA vector

global FID

diff=1E5*(1/SYS.REA.targetKeff-1/SYS.KEFF.EQL0D(end)); % reactivity in pcm

if SYS.REA.allowRemoval&&SYS.REA.allowAddition
  criterion=abs(diff)>abs(SYS.REA.tol);
elseif SYS.REA.allowRemoval&&~SYS.REA.allowAddition
  criterion=diff>abs(SYS.REA.tol);
elseif SYS.REA.allowAddition&&~SYS.REA.allowRemoval
  criterion=diff<-abs(SYS.REA.tol);
end

if criterion %activate reactivity control if above tolerance
  saveTarget=MAT(SYS.IDX.REA.target).N(:,end);
  if ~isempty(SYS.REA.feedMat)
    saveFeed=MAT(SYS.IDX.REA.feed).N(:,end);
  end
  fprintf(FID.log,'%s\n',['** REACT ** k-eff: ' num2str(SYS.KEFF.EQL0D(end)) ' > tol. limit, searching...']);

  switch SYS.REA.mode
    case {'replace','addMass'}
      j=0;
      changeUp=0;
      if diff(1)>0
        maxBound=0; minBound=-1; changeUp(end+1)=-0.99;replReset=false;
        if isempty(SYS.REA.replFraction)&&strcmp(SYS.REA.mode,'replace')
          replReset=true; SYS.REA.replFraction=MAT(SYS.IDX.REA.feed).mFrac(SYS.IDX.REA.feedNucRepl);
        end
      elseif diff(1)<0
        maxBound=1; minBound=0; changeUp(end+1)=0.99; replReset=false;
        if isempty(SYS.REA.replFraction)&&strcmp(SYS.REA.mode,'replace')
          replReset=true; SYS.REA.replFraction=MAT(SYS.IDX.REA.target).mFrac(SYS.IDX.REA.targetNucRepl);
        end
      end
      upReset=false;
      if isempty(SYS.REA.upFraction)  % mass fractions given by feed material
        upReset=true; SYS.REA.upFraction=MAT(SYS.IDX.REA.feed).mFrac(SYS.IDX.REA.feedNucUp);
      end
      downReset=false;
      if isempty(SYS.REA.downFraction)
        downReset=true; SYS.REA.downFraction=MAT(SYS.IDX.REA.target).mFrac(SYS.IDX.REA.targetNucDo);
      end
      while j<SYS.REA.maxIter&&abs(diff(end))>SYS.REA.tol
        j=j+1;
        if j==1
          MAT(SYS.IDX.REA.target).N(:,end+1)=saveTarget;
          if ~isempty(SYS.REA.feedMat)
            MAT(SYS.IDX.REA.feed).N(:,end+1)=saveFeed;
          end
        elseif j>1
          MAT(SYS.IDX.REA.target).N(:,end)=saveTarget;
          if ~isempty(SYS.REA.feedMat)
            MAT(SYS.IDX.REA.feed).N(:,end)=saveFeed;
          end
          if j==SYS.REA.maxIter
            changeUp(end+1)=changeUp(max(find(abs(diff)==min(abs(diff)))));
          else
            changeUp(end+1)=max([min([(changeUp(end-1)*diff(end)-changeUp(end)*diff(end-1))/...
              (diff(end)-diff(end-1)) maxBound]),minBound]);
          end
        end
        NChange=0;
        if diff(1)>0
          NChange=changeUp(end)*MAT(SYS.IDX.REA.target).N(SYS.IDX.REA.targetNucDo,end);
          MAT(SYS.IDX.REA.target).N(SYS.IDX.REA.targetNucDo,end)=...
            MAT(SYS.IDX.REA.target).N(SYS.IDX.REA.targetNucDo,end)+NChange;
          if strcmp(SYS.REA.mode,'replace')
            NChangeRepl=SYS.REA.replFraction.*sum(NChange.*MAT(SYS.IDX.REA.target).atomicMass(SYS.IDX.REA.targetNucDo))...
              ./MAT(SYS.IDX.REA.target).avMass(SYS.IDX.REA.targetNucRepl);
            MAT(SYS.IDX.REA.target).N(SYS.IDX.REA.targetNucRepl,end)=...
              MAT(SYS.IDX.REA.target).N(SYS.IDX.REA.targetNucRepl,end)-NChangeRepl;
          end
        elseif diff(1)<0
          NChange=changeUp(end)*MAT(SYS.IDX.REA.feed).N(SYS.IDX.REA.feedNucUp,end);
          MAT(SYS.IDX.REA.target).N(SYS.IDX.REA.targetNucUp,end)=...
            MAT(SYS.IDX.REA.target).N(SYS.IDX.REA.targetNucUp,end)+NChange;
          if strcmp(SYS.REA.mode,'replace')
            NChangeRepl=SYS.REA.replFraction.*sum(NChange.*MAT(SYS.IDX.REA.feed).atomicMass(SYS.IDX.REA.feedNucUp))...
              ./MAT(SYS.IDX.REA.feed).avMass(SYS.IDX.REA.feedNucRepl);
            MAT(SYS.IDX.REA.target).N(SYS.IDX.REA.targetNucRepl,end)=...
              MAT(SYS.IDX.REA.target).N(SYS.IDX.REA.targetNucRepl,end)-NChangeRepl;
          end
        end
        if ~isempty(SYS.REA.feedMat)
          if diff(1)>0
            MAT(SYS.IDX.REA.feed).N(SYS.IDX.REA.feedNucDo,end)=...
              MAT(SYS.IDX.REA.feed).N(SYS.IDX.REA.feedNucDo,end)-NChange;
          elseif diff(1)<0
            MAT(SYS.IDX.REA.feed).N(SYS.IDX.REA.feedNucUp,end)=...
              MAT(SYS.IDX.REA.feed).N(SYS.IDX.REA.feedNucUp,end)-NChange;
          end
          if strcmp(SYS.REA.mode,'replace')
            MAT(SYS.IDX.REA.feed).N(SYS.IDX.REA.feedNucRepl,end)=...
              MAT(SYS.IDX.REA.feed).N(SYS.IDX.REA.feedNucRepl,end)+NChangeRepl;
          end
        end
        [keff,~]=computeK(MAT,SYS);
        diff(end+1)=1e5*(1/keff-1/SYS.REA.targetKeff);
        if SYS.verboseMode
          fprintf(FID.log,'%s\n',['** REACT ** Current k-eff: ' num2str(SYS.KEFF.EQL0D(end))]);
        end
      end
      if diff(1)>0
        NChange=NChange.*MAT(SYS.IDX.REA.target).atomicMass(SYS.IDX.REA.targetNucDo)*1E24;
      elseif diff(1)<0
        NChange=NChange.*MAT(SYS.IDX.REA.target).atomicMass(SYS.IDX.REA.targetNucUp)*1E24;
      end
      if strcmp(SYS.REA.mode,'addMass')
        elemIdx=unique([SYS.IDX.REA.targetNucUp;SYS.IDX.REA.targetNucDo]);
        NChangeWrite=zeros(size(elemIdx));
        if diff(1)>0
          NChangeWrite(ismember(elemIdx,SYS.IDX.REA.targetNucDo))=NChange;
        elseif diff(1)<0
          NChangeWrite(ismember(elemIdx,SYS.IDX.REA.targetNucUp))=NChange;
        end
      elseif strcmp(SYS.REA.mode,'replace')
        elemIdx=unique([SYS.IDX.REA.targetNucRepl;SYS.IDX.REA.targetNucUp;SYS.IDX.REA.targetNucDo]);
        NChangeWrite=zeros(size(elemIdx));
        if diff(1)>0
          NChangeWrite(ismember(elemIdx,SYS.IDX.REA.targetNucDo))=NChange;
        elseif diff(1)<0
          NChangeWrite(ismember(elemIdx,SYS.IDX.REA.targetNucUp))=NChange;
        end
        NChangeRepl=-NChangeRepl.*MAT(SYS.IDX.REA.target).atomicMass(SYS.IDX.REA.feedNucRepl)*1E24;
        NChangeWrite(ismember(elemIdx,SYS.IDX.REA.targetNucRepl))=NChangeRepl;
      end
      if ~isempty(SYS.IDX.REA.feed)
        SYS.reactAdditions=[NChangeWrite' -NChangeWrite']/1000;
      else
        SYS.reactAdditions=NChangeWrite'/1000;
      end
      if upReset
        SYS.REA.upFraction=[];
      end
      if downReset
        SYS.REA.downFraction=[];
      end
      if replReset
        SYS.REA.replFraction=[];
      end
      if abs(changeUp)==1
        fprintf(FID.log,'%s\n','** REACT ** Max. composition change reached!');
      end
      if ~SYS.PCC.active||SYS.PCC.corrector
        fprintf(FID.log,'%s\n',['** REACT ** Time: ' num2str(SYS.nowTime(end))...
          ' EFPD Change: ' num2str(sum(NChange/1000.0),'%8.6G') ' kg.']);
        fprintf(FID.react,'%-7.3d%-6.3d%-9G',SYS.ouCntr,SYS.inCntr,SYS.nowTime(end));
        fprintf(FID.react,[repmat('%-13.6G',1,numel(SYS.reactAdditions)) '\n'],SYS.reactAdditions);
      end
    case 'addVolume'
      j=0;
      saveVol=MAT(SYS.IDX.REA.target).volume;
      addVolume=[];
      changeVol=0;
      if diff(1)<0
        maxBound=0; minBound=-1; changeVol(end+1)=-0.05;
      elseif diff(1)>0
        maxBound=1; minBound=0; changeVol(end+1)=0.05;
      end
      while j<SYS.REA.maxIter&&abs(diff(end))>SYS.REA.tol
        j=j+1;
        if j==1
          MAT(SYS.IDX.REA.target).N(:,end+1)=saveTarget;
          if ~isempty(SYS.REA.feedMat)
            MAT(SYS.IDX.REA.feed).N(:,end+1)=saveFeed;
          end
        elseif j>1
          MAT(SYS.IDX.REA.target).N(:,end)=saveTarget;
          MAT(SYS.IDX.REA.target).volume=saveVol;
          if j==SYS.REA.maxIter
            changeVol(end+1)=changeVol(max(find(abs(diff)==min(abs(diff)))));
          else
            changeVol(end+1)=max([min([(changeVol(end-1)*diff(end)-changeVol(end)*diff(end-1))/...
              (diff(end)-diff(end-1)) maxBound]), minBound]);
          end
        end
        addVolume(end+1)=changeVol(end)*saveVol;
        MAT(SYS.IDX.REA.target).N(SYS.IDX.targetNuc,end)=...
          MAT(SYS.IDX.REA.target).N(SYS.IDX.targetNuc,end)+addVolume(end)*MAT(SYS.IDX.REA.feed).atDens(SYS.IDX.feedNuc);
        MAT(SYS.IDX.REA.target).volume=MAT(SYS.IDX.REA.target).volume+addVolume(end);
        [keff,~]=computeK(MAT,SYS);
        diff(end+1)=1e5*(1/keff-1/SYS.REA.targetKeff);
      end
      MAT(SYS.IDX.REA.feed).N(SYS.IDX.feedNuc,end)=...
        MAT(SYS.IDX.REA.feed).N(SYS.IDX.feedNuc,end)-addVolume(end)*MAT(SYS.IDX.REA.feed).atDens(SYS.IDX.feedNuc);
      if SYS.verboseMode
        fprintf(FID.log,'%s\n',['** REACT ** Current k-eff: ' num2str(SYS.KEFF.EQL0D(end))]);
      end
      if ~SYS.PCC.active||SYS.PCC.corrector
        MAT(SYS.IDX.REA.feed).volume=MAT(SYS.IDX.REA.feed).volume-addVolume(end);
        fprintf(FID.volume,'%-7.3d %-6.3d %-9G %#-13.6G\n',SYS.ouCntr,SYS.inCntr,SYS.nowTime(end),addVolume(end));
      end
  end
  if j>=SYS.REA.maxIter
    fprintf(FID.log,'%s\n','** REACT ** Max. # iterations reached!');
  else
    fprintf(FID.log,'%s\n',['** REACT ** k-eff: ' num2str(SYS.KEFF.EQL0D(end)) ', stopping... ']);
  end
else
  fprintf(FID.log,'%s\n',['** REACT ** k-eff: ' num2str(SYS.KEFF.EQL0D(end)) ' < tol. limit, no search.']);
end
return
end
