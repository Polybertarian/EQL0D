function [MAT,SYS] = loadSerpentData(MAT,SYS)
    % [MAT,SYS] = LOADSERPENTDATA(MAT,SYS) Loads data from the Serpent outputs
    %in the current folder, such as:
    %   - burn-up matrices from the depmtx_<material>0.m files
    %   - integral fluxes from the <case>_det0.m file
    %   - nubar, leakage probability, k-inf and k-eff as computed by Serpent
    %       from the <case>_res.m file
    %   - nuclide reaction rates from the <case>_arr0.m file and processes them
    %       into the materials vector MAT and system data structure SYS

    %%% Check for files
    exst=[exist([SYS.Casename '_arr0.m'],'file')==2,...
    exist([SYS.Casename  '_res.m'],'file')==2,...
    exist([SYS.Casename '_det0.m'],'file')==2];
    for i=SYS.IDX.MAT.burn
        exst(end+1)=(exist(['depmtx_' MAT(i).name '0.m'],'file')==2);
    end
    if ~all(exst)
        error('Error: Serpent outputs not found!')
    end

    %%% Remove oldest step data
    SYS.RR.notInMat(1)=[]; SYS.RR.inMat(1)=[];
    SYS.RR.NU(1,:)=[]; SYS.RR.LEAK(1)=[];
    %SYS.RR.intCapt(1,:)=[]; SYS.RR.intCapt(3,1:numel(MAT))=zeros(1,numel(MAT));
    %SYS.RR.intFiss(1,:)=[]; SYS.RR.intFiss(3,1:numel(MAT))=zeros(1,numel(MAT));
    %SYS.RR.intProd(1,:)=[]; SYS.RR.intProd(3,1:numel(MAT))=zeros(1,numel(MAT));
    %SYS.RR.LEAK(3)=0;
    SYS.RR.NU(3,:)=zeros(1,numel(SYS.IDX.MAT.inFlux));
    SYS.RR.devFiss(1,:)=[];
    SYS.RR.devCapt(1,:)=[];
    SYS.RR.devFiss(3,:)=zeros(1,numel(SYS.IDX.MAT.inFlux));
    SYS.RR.devCapt(3,:)=zeros(1,numel(SYS.IDX.MAT.inFlux));

    %%% Get integral fluxes
    run([SYS.Casename '_det0.m']);
    for i=SYS.IDX.MAT.inFlux
        MAT(i).intFlux=DETintFlux(i,11);
        SYS.RR.intCapt(i)=DETintCapt(i,11);
        SYS.RR.intFiss(i)=DETintFiss(i,11);
        SYS.RR.intProd(i)=DETintProd(i,11);
        if SYS.RR.intFiss(i)~=0
            SYS.RR.NU(3,i)=SYS.RR.intProd(i)/SYS.RR.intFiss(i);
        end
    end
    SYS.intFlux=sum([MAT(SYS.IDX.MAT.inFlux).intFlux]);
    SYS.tgtFissRate=sum(SYS.RR.intFiss);

    %%% Read data for burnable materials
    for i=SYS.IDX.MAT.burn
        MAT(i)=MAT(i).loadBurnMatrix;
    end

    %%% Get one-group XS/Nubar/Serpent k-eff/inf
    run([SYS.Casename '_res.m']);
    idxUni=find(all(ismember(GC_UNIVERSE_NAME,'0'),2));
    idxUni=idxUni(1);
    SYS.RR.LEAK(3) =ABS_KINF(idxUni,1)/ABS_KEFF(idxUni,1);
    SYS.KEFF.Serpent(end+1)=ABS_KEFF(idxUni,1);
    SYS.KINF.Serpent(end+1)=ABS_KINF(idxUni,1);

    %%% Get micro RR
    run([SYS.Casename '_arr0.m']);
    rr(:,[1 4 5 7])=[]; %remove useless data
    rr(rr(:,3)==0,:)=[]; % remove zero rates
    ZAI=unique(rr(:,1));
    rr_sorted=zeros(length(ZAI),4);
    for i=1:length(ZAI)
        rr_sorted(i,1)=sum(rr(rr(:,1)==ZAI(i)&ismember(rr(:,2),[102:1:117 600:1:849]),3));%capture (n,0n)
        rr_sorted(i,2)=sum(rr(rr(:,1)==ZAI(i)&ismember(rr(:,2),[18:1:21 38]),3));%(n,fission)
        rr_sorted(i,3)=sum(rr(rr(:,1)==ZAI(i)&ismember(rr(:,2),[11 16 24 30 41 875:1:891]),3));%(n,2n)
        rr_sorted(i,4)=sum(rr(rr(:,1)==ZAI(i)&ismember(rr(:,2),[17 25 42]),3));%(n,3n)
    end
    RROidx=~ismember(ZAI,MAT(1).ZAI); %not in Mat
    SYS.RR.notInMat{3}=struct('ZAI',num2cell([-1;ZAI(RROidx)]),'capt',num2cell([0;rr_sorted(RROidx,1)]),...
    'fiss',num2cell([0;rr_sorted(RROidx,2)]),'n2n',num2cell([0;rr_sorted(RROidx,3)]),...
    'n3n',num2cell([0;rr_sorted(RROidx,4)]));
    rr_sorted(RROidx,:)=[]; ZAI(RROidx)=[];

    RR=zeros(length(MAT(1).ZAI),4);
    ZAIidx=ismember(MAT(1).ZAI,ZAI); %in Mat
    RR(ZAIidx,:)=rr_sorted;
    NPhi=[];
    for i=SYS.IDX.MAT.inFlux
        NPhi(:,i)=MAT(i).intFlux*MAT(i).atDens;
    end
    RR=RR./repmat(sum(NPhi,2),1,4); % create cross-sections using integral flux and nuclides densities
    RR(isnan(RR)|isinf(RR))=0.0;
    SYS.RR.inMat{3}=struct('capt',RR(:,1),'fiss',RR(:,2),'n2n',RR(:,3),'n3n',RR(:,4));
    for i=SYS.IDX.MAT.inFlux
        SYS.RR.devFiss(3,i)=SYS.RR.intFiss(i)/sum(NPhi(:,i).*RR(:,2));
        SYS.RR.devCapt(3,i)=SYS.RR.intCapt(i)/sum(NPhi(:,i).*RR(:,1));
    end
    SYS.RR.devFiss(3,isnan(SYS.RR.devFiss(3,:)))=0.0;
    SYS.RR.devCapt(3,isnan(SYS.RR.devCapt(3,:)))=0.0;

    MAT = updateRates(MAT,SYS);

    return
end

%sigma_1*phi1*N1+sigma_2*phi2*N2 = ftot
%sigma_av = ftot/(phi1*N1+phi2*N2)
%sigma_1 = F1/(phi1*N1)
