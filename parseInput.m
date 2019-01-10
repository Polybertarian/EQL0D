function [MAT,OPT,REP,SYS] = parseInput(MAT,OPT,REP,SYS)
    %[MAT,OPT,REP,SYS] = PARSEINPUT(MAT,OPT,REP,SYS) parses the input data

    if SYS.debugMode  %force keeping files in debug mode
        OPT.keepFiles=true;
    end

    if OPT.PCC
        fprintf('%s\n','****  PCC  **** Predictor-Corrector is ON.');
        OPT.renormalize=false;
    else
        fprintf('%s\n','****  PCC  **** Predictor-Corrector is OFF.');
    end

    % Adapt vectors (Cycle length, etc.)
    if length(OPT.cycleLength)<OPT.nCycles  % repeat last value in cycleLength vector to match nCycles
        OPT.cycleLength=[OPT.cycleLength repmat(OPT.cycleLength(end),1,OPT.nCycles-length(OPT.cycleLength))];
    end
    if length(OPT.nSteps)<OPT.nCycles   % repeat last value in nSteps vector to match nCycles
        OPT.nSteps=[OPT.nSteps repmat(OPT.nSteps(end),1,OPT.nCycles-length(OPT.nSteps))];
    end

    %% Reactivity control
    if OPT.reactControl
        fprintf('%s\n','**** REACT **** Reactivity control is ON.');
        SYS.REA=OPT.REA;
        SYS.reactControl=true;
        if ~iscolumn(OPT.REA.upFraction)
            SYS.REA.upFraction=OPT.REA.upFraction';
        end
        if ~iscolumn(OPT.REA.downFraction)
            SYS.REA.downFraction=OPT.REA.downFraction';
        end
        if ~iscolumn(OPT.REA.replFraction)
            SYS.REA.replFraction=OPT.REA.replFraction';
        end

        SYS.IDX.REA.targetMat=find(strcmp({MAT.name},OPT.REA.targetMat));
        SYS.IDX.REA.upNucl=MAT(SYS.IDX.REA.targetMat).find(OPT.REA.upNuclides);
        SYS.IDX.REA.downNucl=MAT(SYS.IDX.REA.targetMat).find(OPT.REA.downNuclides);
        SYS.IDX.REA.replNuc=MAT(SYS.IDX.REA.targetMat).find(OPT.REA.replNuclides);
        if isempty(SYS.IDX.REA.targetMat)
            error('EQL0D:ReactControlTargetMatNotFound',['Error: Reactivity control target material ''' OPT.REA.targetMat ''' not found!']);
        end
        if ~isempty(OPT.REA.feedMat)
            SYS.IDX.REA.feedMat=find(strcmp({MAT.name},OPT.REA.feedMat));
            if isempty(SYS.IDX.REA.feedMat)
                error('EQL0D:ReactControlSourceMatNotFound',['Error: Reactivity control source material ''' OPT.REA.feedMat ''' not found!']);
            end
        end
        if any(OPT.REA.upNuclides<999)
            SYS.REA.upFraction=[];
        end
        if any(OPT.REA.downNuclides<999)
            SYS.REA.downFraction=[];
        end
        SYS.reactAdditions=[];
        SYS.REA.upNuclides=MAT(SYS.IDX.REA.targetMat).ZAI(SYS.IDX.REA.upNucl);
        SYS.REA.downNuclides=MAT(SYS.IDX.REA.targetMat).ZAI(SYS.IDX.REA.downNucl);
    else
        fprintf('%s\n','**** REACT **** Reactivity control is OFF.');
        SYS.IDX.REA=[];
        SYS.REA=[];
        SYS.reactControl=false;
    end

    %% Reprocessing
    if ~isempty(REP)
        matsNeeded=unique({REP.srcMat REP.dstMat}); % materials declared as source/destination for streams
        matsExisting={'' MAT.name}; % materials declared in input
        for j=find(~ismember(matsNeeded,matsExisting))
            MAT(end+1)=Mat(matsNeeded{j},0,0,1e6,0,[],[]); % create dummy materials
        end

        REP=REP([find([REP.isCont]&~[REP.isKeep]) find([REP.isCont]&[REP.isKeep])... % Look for continuous/batch reprocessing streams
        find([REP.isBatch]&~[REP.isKeep]) find([REP.isBatch]&[REP.isKeep])]);

        SYS.IDX.REP.cont =find([REP.isCont]);
        SYS.IDX.REP.batch=find([REP.isBatch]);
        SYS.IDX.REP.keep =find([REP.isKeep]);
        SYS.IDX.REP.contKeep=find([REP.isKeep]&[REP.isCont]);
        for i=1:length(REP) % Identify source and destination material indexes
            REP(i).srcMatIdx=REP(i).findSrc({MAT.name});
            REP(i).dstMatIdx=REP(i).findDst({MAT.name});
            REP(i)=REP(i).adaptElements(MAT(1).ZAI);
            if ~isempty(REP(i).srcMatIdx)
                REP(i).srcNucIdx=ismember(MAT(REP(i).srcMatIdx).ZAI,REP(i).elements);
                if REP(i).isKeep
                    REP(i).srcNucMass=MAT(REP(i).srcMatIdx).atomicMass(REP(i).srcNucIdx);
                end
            end
            if ~isempty(REP(i).dstMatIdx)
                REP(i).dstNucIdx=ismember(MAT(REP(i).dstMatIdx).ZAI,REP(i).elements);
                if REP(i).isKeep
                    REP(i).dstNucMass=MAT(REP(i).dstMatIdx).atomicMass(REP(i).dstNucIdx);
                end
            end
        end
        for i=SYS.IDX.REP.batch % find nuclide positions affected by batch streams
            if isempty(REP(i).srcMatIdx)&&~isempty(REP(i).dstMatIdx)
                REP(i).srcNucIdx=false(size(MAT(REP(i).dstMatIdx).ZAI));
                REP(i).dstNucIdx=ismember(MAT(REP(i).dstMatIdx).ZAI,REP(i).elements);
            elseif isempty(REP(i).dstMatIdx)&&~isempty(REP(i).srcMatIdx)
                REP(i).dstNucIdx=false(size(MAT(REP(i).srcMatIdx).ZAI));
                REP(i).srcNucIdx=ismember(MAT(REP(i).srcMatIdx).ZAI,REP(i).elements);
            elseif ~isempty(REP(i).dstMatIdx)&&~isempty(REP(i).srcMatIdx)
                REP(i).srcNucIdx=ismember(MAT(REP(i).srcMatIdx).ZAI,REP(i).elements);
                REP(i).dstNucIdx=ismember(MAT(REP(i).dstMatIdx).ZAI,MAT(REP(i).srcMatIdx).ZAI(REP(i).srcNucIdx));
            else
                error('Error: Source and/or Destination not recognized while processing')
            end
        end
        for j=SYS.IDX.REP.batch % Convert to mass fraction for batch streams
            if ismember(REP(j).mode,{'keepTotM','keepAM','keepAFPM'})
                REP(j).share=REP(j).share.*MAT(REP(j).dstMatIdx).atomicMass(REP(j).dstNucIdx);
                REP(j).share=REP(j).share/sum(REP(j).share);
            end
        end
        SYS.IDX.MAT.cont=find(ismember({MAT.name},unique({REP([REP.isCont]).srcMat REP([REP.isCont]).dstMat}))); % indexes refering to materials in flux
    else %no reprocessing
        SYS.IDX.REP.cont=[];
        SYS.IDX.REP.batch=[];
        SYS.IDX.MAT.cont=[];
    end
    SYS.IDX.MAT.burn=find([MAT.isBurned]);        % Index of materials burned in flux
    SYS.IDX.MAT.inFlux=find([MAT.isInFlux]);      % index of materials in flux but not burned
    SYS.IDX.MAT.decay=find(~[MAT.isInFlux]&[MAT.isCont]); % Index of decaying materials (out of flux)

    SYS.IDX.MAT.contStreams={};
    for i=SYS.IDX.MAT.cont
        SYS.IDX.MAT.contStreams{end+1}=[];
        for j=SYS.IDX.REP.cont
            if strcmp(REP(j).srcMat,MAT(i).name)||strcmp(REP(j).dstMat,MAT(i).name)
                SYS.IDX.MAT.contStreams{end}(end+1)=j;
            end
        end
    end

    SYS.IDX.REP.matGroups={}; % material groups
    i=unique([SYS.IDX.MAT.burn SYS.IDX.MAT.cont]);
    while ~isempty(i)
        SYS.IDX.REP.matGroups{end+1}=i(1);
        for j=SYS.IDX.REP.cont
            if strcmp(REP(j).srcMat,MAT(i(1)).name)
                SYS.IDX.REP.matGroups{end}=[SYS.IDX.REP.matGroups{end} REP(j).dstMatIdx];
            elseif strcmp(REP(j).dstMat,MAT(i(1)).name)
                SYS.IDX.REP.matGroups{end}=[SYS.IDX.REP.matGroups{end} REP(j).srcMatIdx];
            end
        end
        SYS.IDX.REP.matGroups{end}=sort(unique(SYS.IDX.REP.matGroups{end}));
        i(ismember(i,SYS.IDX.REP.matGroups{end}))=[];
    end
    SYS.IDX.REP.matGroups=[SYS.IDX.REP.matGroups num2cell(SYS.IDX.MAT.decay(~ismember(SYS.IDX.MAT.decay,SYS.IDX.MAT.cont)))];

    for i=1:length(REP) % streams for each material
        switch REP(i).isCont
        case true
            if ~isempty(REP(i).srcMatIdx)
                MAT(REP(i).srcMatIdx).streams.cont(end+1)=i;
            end
            if ~isempty(REP(i).dstMatIdx)
                MAT(REP(i).dstMatIdx).streams.cont(end+1)=i;
            end
        case false
            if ~isempty(REP(i).srcMatIdx)
                MAT(REP(i).srcMatIdx).streams.batch(end+1)=i;
            end
            if ~isempty(REP(i).dstMatIdx)
                MAT(REP(i).dstMatIdx).streams.batch(end+1)=i;
            end
        end
    end

    %% Redox control material indexes
    if OPT.redoxControl
        fprintf('%s\n','**** REDOX **** Redox control is ON.');
        SYS.IDX.redoxMat=find(strcmp({MAT.name},OPT.REDOX.materials));
        SYS.IDX.redoxHalide={MAT(SYS.IDX.redoxMat).mainHalide};
        if isfield(OPT.REDOX,'changeElement')
            for i=SYS.IDX.redoxMat
                for j=1:length(OPT.REDOX.changeElement(:,1))
                    MAT(i).oxState(MAT(i).find(OPT.REDOX.changeElement(j,1)))=OPT.REDOX.changeElement(j,2);
                end
            end
        end
        if ~isempty(OPT.REDOX.replaceMode)&&~isempty(OPT.REDOX.replaceWith)
            for i=SYS.IDX.redoxMat
                SYS.IDX.redoxNuc{i}=MAT(1).find(OPT.REDOX.replaceWith);
            end
        end
    else
        fprintf('%s\n','**** REDOX **** Redox control is OFF.');
    end

    for i=SYS.IDX.MAT.burn     % Adapt things for materials outside of the flux
        [SYS.MTX.decay{[1 2],i}]=deal([]);
        [SYS.MTX.burn{[1 2],i}]=deal(spalloc(length(SYS.MTX.decay{1,i}),length(SYS.MTX.decay{1,i}),0));
        [SYS.IDX.matZAI{[1 2],i}]=deal([]);
        [SYS.IDX.burnZAI{[1 2],i}]=deal([]);
        [SYS.burnZAI{[1 2],i}]=deal(MAT(i).ZAI);
    end
    for i=SYS.IDX.REP.cont
        [SYS.MTX.rep{[1 2],i}]=deal([]);
    end
    [SYS.MTX.total{[1 2]}]=deal([]);

end
