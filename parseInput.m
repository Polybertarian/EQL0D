function [MAT,OPT,REP,SYS] = parseInput(MAT,OPT,REP,SYS)
%[MAT,OPT,REP,SYS] = PARSEINPUT(MAT,OPT,REP,SYS) parses the input data
% Adapt vectors (Cycle length, etc.)
if(length(OPT.cycleLength)<OPT.nCycles) % repeat last value in cycleLength vector to match nCycles
    OPT.cycleLength=[OPT.cycleLength repmat(OPT.cycleLength(end),1,OPT.nCycles-length(OPT.cycleLength))];
end
if(length(OPT.nSteps)<OPT.nCycles)  % repeat last value in nSteps vector to match nCycles
    OPT.nSteps=[OPT.nSteps repmat(OPT.nSteps(end),1,OPT.nCycles-length(OPT.nSteps))];
end

%% Reactivity control
if(OPT.reactControl)
    if ~iscolumn(OPT.REA.upFraction)
        OPT.REA.upFraction=OPT.REA.upFraction';
    end
    if ~iscolumn(OPT.REA.downFraction)
        OPT.REA.downFraction=OPT.REA.downFraction';
    end
    if ~iscolumn(OPT.REA.replFraction)
        OPT.REA.replFraction=OPT.REA.replFraction';
    end
    
    SYS.IDX.targetMat=find(strcmp({MAT.name},OPT.REA.targetMat));
    if(~isempty(OPT.REA.feedMat))
        SYS.IDX.feedMat=find(strcmp({MAT.name},OPT.REA.feedMat));
        if(strcmp(OPT.REA.mode,'addVolume'))
            SYS.IDX.feedNuc=MAT(SYS.IDX.feedMat).find(OPT.REA.upNuclides);
        else
            SYS.IDX.feedNucUp=MAT(SYS.IDX.feedMat).find(OPT.REA.upNuclides);
            SYS.IDX.feedNucDo=MAT(SYS.IDX.feedMat).find(OPT.REA.downNuclides);
            if(ismember(OPT.REA.mode,{'replace'}))
                SYS.IDX.feedNucRepl=MAT(SYS.IDX.feedMat).find(OPT.REA.replNuclides);
            end
        end
    end
    if(strcmp(OPT.REA.mode,'addVolume'))
        SYS.IDX.targetNuc=MAT(SYS.IDX.targetMat).find(OPT.REA.upNuclides);
    else
        SYS.IDX.targetNucUp=MAT(SYS.IDX.targetMat).find(OPT.REA.upNuclides);
        SYS.IDX.targetNucDo=MAT(SYS.IDX.targetMat).find(OPT.REA.downNuclides);
        if(ismember(OPT.REA.mode,{'replace'}))
            SYS.IDX.targetNucRepl=MAT(SYS.IDX.targetMat).find(OPT.REA.replNuclides);
            if(any(OPT.REA.replNuclides<999))
                OPT.REA.replFraction=[];
            end
        end
    end
    if(any(OPT.REA.upNuclides<999))
        OPT.REA.upFraction=[];
    end
    if(any(OPT.REA.downNuclides<999))
        OPT.REA.downFraction=[];
    end
    SYS.reactAdditions=[];
end

%% Reprocessing
if(~isempty(REP))
    %%% Create 'dummy' materials for streams
    matsNeeded=unique({REP.srcMat REP.dstMat}); % materials declared as source/destination for streams
    matsExisting={'' MAT.name}; % materials declared in input
    for j=find(~ismember(matsNeeded,matsExisting))
        MAT(end+1)=Mat(matsNeeded{j},0,0,1e6,0,[],[]); % create dummy materials
    end
    
    %%% Look for continuous/batch reprocessing streams
    REP=REP([find([REP.isCont]&~[REP.isKeep])...
        find([REP.isCont]&[REP.isKeep])...
        find([REP.isBatch]&~[REP.isKeep])...
        find([REP.isBatch]&[REP.isKeep])]);
    
    SYS.IDX.contStr=find([REP.isCont]);
    SYS.IDX.batchStr=find([REP.isBatch]);
    
    %%% Identify source and destination material indexes
    for i=1:length(REP)
        REP(i).srcMatIdx=REP(i).findSrc({MAT.name});
        REP(i).dstMatIdx=REP(i).findDst({MAT.name});
        REP(i)=REP(i).adaptElements(MAT(1).ZAI);
    end
    
    %%% find nuclide positions affected by batch streams
    for i=SYS.IDX.batchStr
        if(isempty(REP(i).srcMatIdx)&~isempty(REP(i).dstMatIdx))
            REP(i).srcNucIdx=false(size(MAT(REP(i).dstMatIdx).ZAI));
            REP(i).dstNucIdx=ismember(MAT(REP(i).dstMatIdx).ZAI,REP(i).elements);
        elseif(isempty(REP(i).dstMatIdx)&~isempty(REP(i).srcMatIdx))
            REP(i).dstNucIdx=false(size(MAT(REP(i).srcMatIdx).ZAI));
            REP(i).srcNucIdx=ismember(MAT(REP(i).srcMatIdx).ZAI,REP(i).elements);
        elseif(~isempty(REP(i).dstMatIdx)&~isempty(REP(i).srcMatIdx))
            REP(i).srcNucIdx=ismember(MAT(REP(i).srcMatIdx).ZAI,REP(i).elements);
            REP(i).dstNucIdx=ismember(MAT(REP(i).dstMatIdx).ZAI,MAT(REP(i).srcMatIdx).ZAI(REP(i).srcNucIdx));
        else
            error('Error: Source and/or Destination not recognized while processing')
        end
    end
    
    %%% Convert to mass fraction for batch streams
    for j=SYS.IDX.batchStr
        if(ismember(REP(j).mode,{'keepTotM','keepAM','keepAFPM'}))
            REP(j).share=REP(j).share.*MAT(REP(j).dstMatIdx).atomicMass(REP(j).dstNucIdx);
            REP(j).share=REP(j).share/sum(REP(j).share);
        end
    end
    %%% indexes refering to materials in flux
    
    SYS.IDX.contMat=find(ismember({MAT.name},unique({REP([REP.isCont]).srcMat REP([REP.isCont]).dstMat})));
    SYS.IDX.strMat=find([MAT.isStr]);
    SYS.IDX.contStrMat=SYS.IDX.strMat(ismember(SYS.IDX.strMat,SYS.IDX.contMat));
else %no reprocessing
    SYS.IDX.contStr=[];
    SYS.IDX.batchStr=[];
    SYS.IDX.contMat=[];
    SYS.IDX.strMat=[];
    SYS.IDX.contStrMat=[];
end
SYS.IDX.burnMat=find([MAT.isBurned]);
SYS.IDX.fluxMat=find([MAT.isInFlux]);

%% Redox control material indexes
if(OPT.redoxControl)
    SYS.IDX.redoxMat=find(strcmp({MAT.name},OPT.REDOX.materials));
    SYS.IDX.redoxHalide={MAT(SYS.IDX.redoxMat).mainHalide};
    if(~isempty(OPT.REDOX.changeElement))
        for i=SYS.IDX.redoxMat
            for j=1:length(OPT.REDOX.changeElement(:,1))
                MAT(i).oxState(MAT(i).find(OPT.REDOX.changeElement(j,1)))=OPT.REDOX.changeElement(j,2);
            end
        end
    end
    if(~isempty(OPT.REDOX.replaceMode)&&~isempty(OPT.REDOX.replaceWith))
        for i=SYS.IDX.redoxMat
            SYS.IDX.redoxNuc{i}=MAT(1).find(OPT.REDOX.replaceWith);
        end
    end
end

%% Adapt things for materials outside of the flux
for i=SYS.IDX.burnMat
    [SYS.MTX.decay{[1 2],i}]=deal([]);
    [SYS.MTX.burn{[1 2],i}]=deal(spalloc(length(SYS.MTX.decay{1,i}),length(SYS.MTX.decay{1,i}),0));
    [SYS.IDX.matZAI{[1 2],i}]=deal([]);
    [SYS.IDX.burnZAI{[1 2],i}]=deal([]);
    [SYS.burnZAI{[1 2],i}]=deal(MAT(i).ZAI);
end
for i=SYS.IDX.strMat
    b=[];
    for j=1:numel(REP)
        if(REP(j).dstMatIdx==i|REP(j).srcMatIdx==i)
            b=unique(vertcat(b,REP(j).elements));
        end
    end
    j=0;
    while(j<length(b))
        j=j+1;
        b=unique([b;MAT(i).ZAI(SYS.MTX.defaultDecay(:,MAT(i).ZAI==b(j))>0)],'stable');
    end
    b=sort(b);
    b=MAT(1).ZAI;
    A=SYS.MTX.defaultDecay;
    [SYS.IDX.burnZAI{[1 2],i}]=deal(ismember(MAT(i).ZAI,b));
    A(~SYS.IDX.burnZAI{2,i},:)=[];
    A(:,~SYS.IDX.burnZAI{2,i})=[];
    [SYS.burnZAI{[1 2],i}]=deal(b);
    [SYS.MTX.burn{[1 2],i}]=deal(spalloc(numel(b),numel(b),1));
    [SYS.MTX.decay{[1 2],i}]=deal(A);
    [SYS.IDX.matZAI{[1 2],i}]=deal(repmat(i,numel(b),1));
    clearvars b
end
for i=SYS.IDX.contStr
    [SYS.MTX.rep{[1 2],i}]=deal([]);
end
[SYS.MTX.total{[1 2]}]=deal([]);


