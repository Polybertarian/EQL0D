classdef Rep
  %REP Reprocessing stream in EQL0D
  
  
  
  properties
    name
    srcMat
    srcMatIdx=[];
    srcNucIdx=[];
    srcNucMass=[];
    dstMat
    dstMatIdx=[];
    dstNucIdx=[];
    dstNucMass=[];
    elements
    elementsNames
    share
    frequency
    rate
    type
    mode
    isCont=false
    isBatch=false
    isKeep=false;
  end
  properties (Dependent)
    cycleTime
    repMtx
  end
  methods
    function obj = Rep(name,srcMat,dstMat,elements,share,rate,type)
      global DAT
      obj.name=name;
      obj.srcMat=srcMat;
      obj.dstMat=dstMat;
      if(isempty(share))
        obj.share=ones(size(elements));
      elseif(length(share)~=length(elements))
        error('Error: share and elements vector length mismatch.');
      else
        obj.share=share;
      end
      if(~iscolumn(obj.share))
        obj.share=obj.share';
      end
      if(~iscolumn(elements))
        elements=elements';
      end
      obj.elementsNames=strjoin(ZAI2Name(elements)');
      if(any(elements<111))
        zailist=DAT.ZAI0(isProduced(DAT.libraryName,DAT.ZAI0));
        [idx1,idx2]=isElement(elements,zailist);
        obj.elements=zailist(idx1);
        obj.share=obj.share(idx2);
      else
        obj.elements=elements;
        
      end
      switch type
        case {'cont','continuous'}
          obj.type='cont';
          obj.isCont=true;
        case {'batch','batchwise'}
          obj.type='batch';
          obj.isBatch=true;
        otherwise
          error('Error while setting REP: type not recognized!');
      end
      obj.frequency=1;
      if(ischar(rate))
        obj.rate=1;
        if(ismember(rate,{'keep','keepAM','keepAFPM','keepTotM','keepAA','keepAFPM','keepTotA'}))
          obj.mode=rate;
          obj.isKeep=true;
        else
          error(['Reprocessing stream mode ' rate ' of stream ' obj.name ' not recognized!']);
        end
      elseif(isnumeric(rate))
        obj.rate=rate;
        obj.mode='remove';
      end
    end
    function srcIdx = findSrc(obj,matList)
      if(isempty(obj.srcMat))
        srcIdx=[];
      else
        srcIdx=find(ismember(matList,obj.srcMat));
      end
    end
    function dstIdx = findDst(obj,matList)
      if(isempty(obj.dstMat))
        dstIdx=[];
      else
        dstIdx=find(ismember(matList,obj.dstMat));
      end
    end
    function obj = adaptElements(obj,ZAIList)
      if(any(obj.elements<111))
        [idx1,idx2]=isElement(obj.elements,ZAIList);
        obj.elements=ZAIList(idx1);
        obj.share=obj.share(idx2);
      end
    end
    function [dstMAT,srcMAT] = batchProcessing(obj,dstMAT,srcMAT,tStep)
      global FID
      %%% 1)determine isotopics, 2)determine total quantity, 3) make change
      if(isempty(obj.srcMatIdx))
        frac=obj.share./dstMAT.atomicMass(obj.dstNucIdx)/1.0E24;
      else
        if(ismember(obj.mode,{'keepAFPM','keepAM','keepTotM'}))
          frac=srcMAT.mFrac(obj.srcNucIdx)./srcMAT.atomicMass(obj.srcNucIdx)/1.0E24;
        elseif(ismember(obj.mode,{'keepAFPA','keepAA','keepTotA'}))
          frac=srcMAT.aFrac(obj.srcNucIdx);
        else
          frac=srcMAT.aFrac(obj.srcNucIdx);
        end
      end
      if(ismember(obj.mode,{'keepAFPM','keepAM','keepTotM'}))
        fillmode='mass';
      elseif(ismember(obj.mode,{'keepAFPA','keepAA','keepTotA'}))
        fillmode='atomic';
      else
        fillmode='atomic';
      end
      switch obj.mode
        case 'keepTotM' %refill total mass
          mDefect=dstMAT.initTotMass-dstMAT.totMass;
        case 'keepTotA' %refill total nuclides
          mDefect=dstMAT.initTotN-dstMAT.totN;
        case 'keepAFPM' %refill up to initial actinide mass - current FP mass
          mDefect=dstMAT.initTotActMass-dstMAT.totActMass-dstMAT.totFPMass;
        case 'keepAFPA' %refill up to initial actinide mass - current FP mass
          mDefect=dstMAT.initTotActN-dstMAT.totActN-dstMAT.totFPN;
        case 'keepAM'   %refill up to initial actinide mass
          mDefect=dstMAT.initTotActMass-dstMAT.totActMass;
        case 'keepAA'   %refill up to initial actinide nuclides
          mDefect=dstMAT.initTotActN-dstMAT.totActN;
        case 'remove'
          if(obj.rate*tStep<1)
            mDefect=obj.rate*tStep*srcMAT.N(obj.srcNucIdx,end); %not a mass
          else
            mDefect=srcMAT.N(obj.srcNucIdx,end); %not a mass
          end
          frac=obj.share;
      end
      switch fillmode
        case 'atomic'
          fprintf(FID.log,'%s\n',['*** BATCH *** Computed atomic defect: ' num2str(sum(mDefect)*1E24,'%E') ' atoms.']);
        case 'mass'
          fprintf(FID.log,'%s\n',['*** BATCH *** Computed mass defect: ' num2str(sum(mDefect)/1000) ' kg.']);
      end
      if(obj.srcMatIdx~=0)
        srcMAT.N(:,end+1)=srcMAT.N(:,end);
        srcMAT.N(obj.srcNucIdx,end)=srcMAT.N(obj.srcNucIdx,end)-frac.*mDefect;
        changeSrc=srcMAT.N(obj.srcNucIdx,end)-srcMAT.N(obj.srcNucIdx,end-1);
        nucNames={srcMAT.nuclideName{obj.srcNucIdx}};
        srcMatName=srcMAT.name;
      else
        changeSrc=zeros(length(find(obj.dstNucIdx)));
        srcMatName='void';
      end
      if(obj.dstMatIdx~=0)
        dstMAT.N(:,end+1)=dstMAT.N(:,end);
        dstMAT.N(obj.dstNucIdx,end)=dstMAT.N(obj.dstNucIdx,end)+frac.*mDefect;
        changeDst=dstMAT.N(obj.dstNucIdx,end)-dstMAT.N(obj.dstNucIdx,end-1);
        nucNames={dstMAT.nuclideName{obj.dstNucIdx}};
        dstMatName=dstMAT.name;
      else
        changeDst=zeros(length(find(obj.srcNucIdx)));
        dstMatName='void';
      end
      fprintf(FID.log,'%s\n','*** BATCH *** Summary of changes: ');
      fprintf(FID.log,'%8s\t\t%13s\t\t%13s\n','Nuc.',srcMatName,dstMatName);
      for i=1:length(changeDst)
        fprintf(FID.log,'%8s\t\t%+13E\t\t%+13E\n',nucNames{i},1E24*changeSrc(i),1E24*changeDst(i));
      end
    end
    function T = get.cycleTime(obj)
      if(strcmp(obj.mode,'remove'))
        T=(1/obj.rate)/24/3600;
      else
        T=[];
      end
    end
    function [] = disp(obj)
      table({obj.srcMat}',{obj.dstMat}',{obj.elementsNames}',{obj.cycleTime}',{obj.mode}',{obj.type}',...
        'VariableNames',{'Source','Destination','Elements','Cycletime','Mode','Type'},...
        'RowNames',{obj.name}')
    end
    function usesMat=usesMat(obj,matName)
      usesMat=strcmp(obj.srcMat,matName)|strcmp(obj.dstMat,matName);
    end
    function repMtx=get.repMtx(obj)
      repMtx=spalloc(length(obj.srcNucIdx),length(obj.srcNucIdx),length(obj.elements));
      if(~obj.isKeep)
        if(isempty(obj.srcMatIdx)) %%% void source = find destination elements
          repMtx(sub2ind(size(repMtx),find(obj.dstNucIdx),find(obj.dstNucIdx)))=+obj.share*obj.rate;
        else %%% non-void source & destination = find both
          repMtx(sub2ind(size(repMtx),find(obj.srcNucIdx),find(obj.srcNucIdx)))=-obj.share*obj.rate;
        end
      end
    end
    function obj=setIdx(obj,srcMatZAI,dstMatZAI)
      if(~isempty(srcMatZAI))
        obj.srcNucIdx=ismember(srcMatZAI,obj.elements);
      end
      if(~isempty(dstMatZAI))
        obj.dstNucIdx=ismember(dstMatZAI,obj.elements);
      end
    end
  end
  
end

