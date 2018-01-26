function [MAT,OPT,REP,SYS] = parseInput(MAT,OPT,REP,SYS)
%[MAT,OPT,REP,SYS] = PARSEINPUT(MAT,OPT,REP,SYS) parses the input data 
% Adapt vectors (Cycle length, etc.)
if(length(OPT.cycleLength)<OPT.nCycles) % repeat last value in cycleLength vector to match nCycles
    OPT.cycleLength=[OPT.cycleLength repmat(OPT.cycleLength(end),1,OPT.nCycles-length(OPT.cycleLength))];
end
if(length(OPT.nSteps)<OPT.nCycles)  % repeat last value in nSteps vector to match nCycles
    OPT.nSteps=[OPT.nSteps repmat(OPT.nSteps(end),1,OPT.nCycles-length(OPT.nSteps))];
end

%%% Reactivity control
if(OPT.reactControl)
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
    SYS.IDX.contStr=find(ismember({REP.type},{'continuous','cont'}));
    SYS.IDX.batchStr=find(ismember({REP.type},{'batchwise','batch'}));
    
    %%% Identify source and destination material indexes
    for i=1:length(REP)
        SYS.IDX.srcMat(i)=REP(i).findSrc({MAT.name});
        SYS.IDX.dstMat(i)=REP(i).findDst({MAT.name});
        REP(i)=REP(i).adaptElements(MAT(1).ZAI);
    end
    
    %%% find nuclide positions affected by batch streams
    for i=SYS.IDX.batchStr
        switch xor(SYS.IDX.srcMat(i),SYS.IDX.dstMat(i))
            case 1 %%% isotopes/shares are given by material composition
                if(SYS.IDX.srcMat(i)==0)
                    SYS.IDX.srcPos{i}=false(size(MAT(SYS.IDX.dstMat(i)).ZAI));
                    SYS.IDX.dstPos{i}=ismember(MAT(SYS.IDX.dstMat(i)).ZAI,REP(i).elements);
                elseif(SYS.IDX.dstMat(i)==0)
                    SYS.IDX.dstPos{i}=false(size(MAT(SYS.IDX.srcMat(i)).ZAI));
                    SYS.IDX.srcPos{i}=ismember(MAT(SYS.IDX.srcMat(i)).ZAI,REP(i).elements);
                end
            case 0 %%% source and destination fixed
                SYS.IDX.srcPos{i}=ismember(MAT(SYS.IDX.srcMat(i)).ZAI,REP(i).elements);
                SYS.IDX.dstPos{i}=ismember(MAT(SYS.IDX.dstMat(i)).ZAI,MAT(SYS.IDX.srcMat(i)).ZAI(SYS.IDX.srcPos{i}));
            otherwise
                error('Error: Source and/or Destination not recognized while processing')
        end
    end
    
    %%% Convert to mass fraction for batch streams
    for j=SYS.IDX.batchStr
        if(ismember(REP(j).mode,{'keepTotM','keepAM','keepAFPM'}))
            REP(j).share=REP(j).share.*MAT(SYS.IDX.dstMat(j)).atomicMass(SYS.IDX.dstPos{j});
            REP(j).share=REP(j).share/sum(REP(j).share);
        end
    end
    
else %no reprocessing
    SYS.IDX.contStr=[];
    SYS.IDX.batchStr=[];
end

%%% indexes refering to materials in flux
SYS.IDX.burnMat=find([MAT.isBurned]);
SYS.IDX.fluxMat=find([MAT.isInFlux]);
SYS.IDX.contMat=find(ismember({MAT.name},unique({REP([REP.isCont]).srcMat REP([REP.isCont]).dstMat})));
SYS.IDX.strMat=find([MAT.isStr]);

%%% Redox control material indexes
if(OPT.redoxControl)
    SYS.IDX.redoxMat=find(strcmp({MAT.name},OPT.REDOX.materials));
    SYS.IDX.redoxHalide={MAT(SYS.IDX.redoxMat).mainHalide};
    if(~isempty(OPT.REDOX.changeElement))
        for i=SYS.IDX.redoxMat
            for j=1:1:length(OPT.REDOX.changeElement(:,1))
                MAT(i).oxState(isElement(OPT.REDOX.changeElement(j,1),MAT(i).ZAI))=...
                    OPT.REDOX.changeElement(j,2);
            end
        end
    end
    if(~isempty(OPT.REDOX.replaceMode)&&~isempty(OPT.REDOX.replaceWith))
        for i=SYS.IDX.redoxMat
            SYS.IDX.redoxNuc{i}=find(ismember(MAT(i).ZAI,OPT.REDOX.replaceWith));
        end
    end
end

%%% Adapt things for materials outside of the flux
for i=SYS.IDX.burnMat
    [SYS.MTX.decay{[1 2],i}]=deal([]);
    [SYS.MTX.burn{[1 2],i}]=deal(spalloc(length(SYS.MTX.decay{1,i}),length(SYS.MTX.decay{1,i}),0));
    [SYS.IDX.matZAI{[1 2],i}]=deal([]);
    [SYS.IDX.burnZAI{[1 2],i}]=deal([]);
    [SYS.burnZAI{[1 2],i}]=deal(MAT(i).ZAI);
end
for i=SYS.IDX.strMat
    b=unique(vertcat(REP(SYS.IDX.dstMat==i|SYS.IDX.srcMat==i).elements));
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
