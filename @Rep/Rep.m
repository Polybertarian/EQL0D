classdef Rep
  %REP Reprocessing stream in EQL0D
  properties
    name
    srcMat
    srcMatIdx=[];
    nucIdx=[];
    srcNucMass=[];
    dstMat
    dstMatIdx=[];
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
      if isempty(share)
        obj.share=ones(size(elements));
      elseif length(share)~=length(elements)
        error('Error: share and elements vector length mismatch.');
      else
        obj.share=share;
      end
      if ~iscolumn(obj.share)
        obj.share=obj.share';
      end
      if ~iscolumn(elements)
        elements=elements';
      end
      obj.elementsNames=strjoin(ZAI2Name(elements)');
      if any(elements<111)
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
      if ischar(rate)
        obj.rate=1;
        if ismember(rate,{'keep','keepAM','keepAFPM','keepTotM','keepAA','keepAFPM','keepTotA'})
          obj.mode=rate;
          obj.isKeep=true;
        else
          error(['Reprocessing stream mode ' rate ' of stream ' obj.name ' not recognized!']);
        end
      elseif isnumeric(rate)
        obj.rate=rate;
        obj.mode='remove';
      end
    end
    function srcIdx = findSrc(obj,matList)
      if isempty(obj.srcMat)
        srcIdx=[];
      else
        srcIdx=find(ismember(matList,obj.srcMat));
      end
    end
    function dstIdx = findDst(obj,matList)
      if isempty(obj.dstMat)
        dstIdx=[];
      else
        dstIdx=find(ismember(matList,obj.dstMat));
      end
    end
    function obj = adaptElements(obj,ZAIList)
      if any(obj.elements<111)
        [idx1,idx2]=isElement(obj.elements,ZAIList);
        obj.elements=ZAIList(idx1);
        obj.share=obj.share(idx2);
      end
    end
    function [dstMAT,srcMAT] = batchProcessing(obj,dstMAT,srcMAT,tStep)
      %%% 1)determine isotopics, 2)determine total quantity, 3) make change
      if isempty(obj.srcMatIdx)
        frac=obj.share./dstMAT.atomicMass(obj.nucIdx)/1.0E24;
      else
        if ismember(obj.mode,{'keepAFPM','keepAM','keepTotM'})
          frac=srcMAT.mFrac(obj.nucIdx)./srcMAT.atomicMass(obj.nucIdx)/1.0E24;
        elseif ismember(obj.mode,{'keepAFPA','keepAA','keepTotA'})
          frac=srcMAT.aFrac(obj.nucIdx);
        else
          frac=srcMAT.aFrac(obj.nucIdx);
        end
      end
      if ismember(obj.mode,{'keepAFPM','keepAM','keepTotM'})
        fillmode='mass';
      elseif ismember(obj.mode,{'keepAFPA','keepAA','keepTotA'})
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
          if obj.rate*tStep<1
            mDefect=obj.rate*tStep*srcMAT.N(obj.nucIdx,end); %not a mass
          else
            mDefect=srcMAT.N(obj.nucIdx,end); %not a mass
          end
          frac=obj.share;
      end
      switch fillmode
        case 'atomic'
          fprintf('%s\n',['  ** BATCH **   Computed atomic defect: ' num2str(sum(mDefect)*1E24,'%E') ' atoms.']);
        case 'mass'
          fprintf('%s\n',['  ** BATCH **   Computed mass defect: ' num2str(sum(mDefect)/1000) ' kg.']);
      end
      if obj.srcMatIdx~=0
        saveSrc=srcMAT.N(obj.nucIdx,end);
      end
      if obj.dstMatIdx~=0
        saveDst=dstMAT.N(obj.nucIdx,end);
      end
      [srcMAT,dstMAT]=transferNuclides(srcMAT,dstMAT,obj.nucIdx,frac.*mDefect);
      if obj.srcMatIdx~=0
        changeSrc=srcMAT.N(obj.nucIdx,end)-saveSrc;
        nucNames={srcMAT.nuclideName{obj.nucIdx}};
        srcMatName=srcMAT.name;
      else
        changeSrc=zeros(length(find(obj.nucIdx)));
        srcMatName='void';
      end
      if obj.dstMatIdx~=0
        changeDst=dstMAT.N(obj.nucIdx,end)-saveDst;
        nucNames={dstMAT.nuclideName{obj.nucIdx}};
        dstMatName=dstMAT.name;
      else
        changeDst=zeros(length(find(obj.nucIdx)));
        dstMatName='void';
      end
      fprintf('%s\n','  ** BATCH **   Summary of changes: ');
      fprintf('%8s\t\t%13s\t\t%13s\n','Nuc.',srcMatName,dstMatName);
      for i=1:length(changeDst)
        fprintf('%8s\t\t%+13E\t\t%+13E\n',nucNames{i},1E24*changeSrc(i),1E24*changeDst(i));
      end
    end
    function T = get.cycleTime(obj)
      if strcmp(obj.mode,'remove')
        T=(1/obj.rate)/24/3600;
      else
        T=[];
      end
    end
    function disp(obj)
      T=table({obj.srcMat}',{obj.dstMat}',{obj.elementsNames}',{obj.cycleTime}',{obj.mode}',{obj.type}',...
        'VariableNames',{'Source','Destination','Elements','Cycletime','Mode','Type'},...
        'RowNames',{obj.name}');
        disp(T)
    end
    function usesMat=usesMat(obj,matName)
      usesMat=strcmp(obj.srcMat,matName)|strcmp(obj.dstMat,matName);
    end
    function repMtx=get.repMtx(obj)
      repMtx=spalloc(length(obj.nucIdx),length(obj.nucIdx),length(obj.elements));
      if ~obj.isKeep
        if isempty(obj.srcMatIdx) % void source = find destination elements
          repMtx(sub2ind(size(repMtx),find(obj.nucIdx),find(obj.nucIdx)))=+obj.share*obj.rate;
        else % non-void source & destination = find both
          repMtx(sub2ind(size(repMtx),find(obj.nucIdx),find(obj.nucIdx)))=-obj.share*obj.rate;
        end
      end
    end
  end
end
