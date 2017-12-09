%%% Adapt vectors (Cycle length, etc.)
if(length(OPT.cycleLength)<OPT.nCycles) % repeat last value in cycleLength vector to match nCycles
    OPT.cycleLength=[OPT.cycleLength repmat(OPT.cycleLength(end),1,OPT.nCycles-length(OPT.cycleLength))];
end
if(length(OPT.nSteps)<OPT.nCycles)  % repeat last value in nSteps vector to match nCycles
    OPT.nSteps=[OPT.nSteps repmat(OPT.nSteps(end),1,OPT.nCycles-length(OPT.nSteps))];
end

%%% Reactivity control
if(OPT.reactControl)
    SYS.IDX.targetMat=find(strcmp({MAT.name},OPT.REA.targetMat));
    if(ismember(OPT.REA.mode,{'replace','addMass'}))
        if(~isempty(OPT.REA.feedMat))
            SYS.IDX.feedMat=find(strcmp({MAT.name},OPT.REA.feedMat));
            SYS.IDX.feedNucUp=MAT(SYS.IDX.feedMat).find(OPT.REA.upNuclides);
 	    if(ismember(OPT.REA.mode,{'replace'}))
        	SYS.IDX.feedNucRepl=MAT(SYS.IDX.feedMat).find(OPT.REA.replNuclides);
 	    end
        end
	if(any(OPT.REA.upNuclides<999))
    	    OPT.REA.upFraction=[];
 	end
        SYS.IDX.targetNucUp=MAT(SYS.IDX.feedMat).find(OPT.REA.upNuclides);
        if(any(OPT.REA.downNuclides<999))
            OPT.REA.downFraction=[];
        end
        SYS.IDX.targetNucDo=MAT(SYS.IDX.targetMat).find(OPT.REA.downNuclides);
        SYS.IDX.feedNucDo=MAT(SYS.IDX.feedMat).find(OPT.REA.downNuclides);
	if(ismember(OPT.REA.mode,{'replace'}))
		SYS.IDX.targetNucRepl=MAT(SYS.IDX.targetMat).find(OPT.REA.replNuclides);
 		if(any(OPT.REA.replNuclides<999))
     			OPT.REA.replFraction=[];
 		end
	end
    elseif(strcmp(OPT.REA.mode,'addVolume'))
        if(~isempty(OPT.REA.feedMat))
            SYS.IDX.feedMat=find(strcmp({MAT.name},OPT.REA.feedMat));
            SYS.IDX.targetNuc=find(ismember(MAT(SYS.IDX.targetMat).ZAI,...
                MAT(SYS.IDX.feedMat).ZAI(MAT(SYS.IDX.feedMat).atDens>0)));
            SYS.IDX.feedNuc=find(MAT(SYS.IDX.feedMat).atDens>0);
        else
            SYS.IDX.targetNuc=MAT(SYS.IDX.targetMat).find(OPT.REA.upNuclides);
        end
    end
    SYS.reactAdditions=[];
end

%% Reprocessing
if(~isempty(REP))
    if(~isfield(REP,'share'))
        tmp=repmat({1},size(REP));
        [REP.share]=deal(tmp{:});
    end
    if(~isfield(REP,'frequency'))
        tmp=repmat({1},size(REP));
        [REP.frequency]=deal(tmp{:});
    end
    
    %%% Check the modes are recognized
    if(any(~ismember({REP.mode},{'remove','keepAFPM','keepAM','keepTotM','keepAFPA','keepAA','keepTotA'})))
        error(['Reprocessing stream mode ' REP(i).mode ' of stream ' REP(i).name ' not recognized!']);
    end
    
    %%% Create 'dummy' materials for streams
    matsNeeded=unique({REP.source REP.destination}); % materials declared as source/destination for streams
    matsExisting={'void' MAT.name}; % materials declared in input
    for j=find(~ismember(matsNeeded,matsExisting))
        MAT(end+1)=Mat(matsNeeded{j},0,0,1e6,0,[],[]); % create dummy materials 
    end
    
    %%% Look for continuous/batch reprocessing streams
    SYS.IDX.contStr=find(ismember({REP.type},{'continuous','cont'}));
    SYS.IDX.batchStr=find(ismember({REP.type},{'batchwise','batch'}));
    
    %%% Identify source and destination material indexes
    for i=1:length(REP)
        if(~iscolumn(REP(i).share))
            REP(i).share=REP(i).share';
        end
        if(~iscolumn(REP(i).elements))
            REP(i).elements=REP(i).elements';
        end
        
        %%% find source/dest material, if any
        SYS.IDX.srcMat(i)=find(ismember({'void' MAT.name},REP(i).source))-1;
        SYS.IDX.dstMat(i)=find(ismember({'void' MAT.name},REP(i).destination))-1;
        
        %%% clarify shares
        if(length(REP(i).share)==1)
            REP(i).share=REP(i).share(1)*ones(size(REP(i).elements));
        elseif(length(REP(i).share)<length(REP(i).elements))
            error('Error: share and elements vector length mismatch.')
        end
        
        %%% transform elemental input into isotopes
        if(any(REP(i).elements<111))
            [idx1,idx2]=isElement(REP(i).elements,MAT(SYS.IDX.srcMat(i)).ZAI);
            REP(i).elements=MAT(SYS.IDX.srcMat(i)).ZAI(idx1);
            REP(i).share=REP(i).share(idx2);
        end
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
SYS.IDX.strMat=find([MAT.isStr]);
SYS.IDX.contMat=find([MAT.isCont]);

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
    if(~isempty(OPT.REDOX.replaceMode)&~isempty(OPT.REDOX.replaceWith))
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
