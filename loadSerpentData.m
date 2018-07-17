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
if(~all(exst))
    error('Error: Serpent outputs not found!')
end

%%% Remove oldest step data
SYS.MTX.decay(1,:)=[]; SYS.MTX.burn(1,:)=[];
SYS.IDX.burnZAI(1,:)=[]; SYS.IDX.matZAI(1,:)=[];
SYS.RR.notInMat(1)=[]; SYS.RR.inMat(1)=[];
SYS.RR.NU(1)=[]; SYS.RR.LEAK(1)=[];
SYS.burnZAI(1,:)=[];

%%% Get integral fluxes
run([SYS.Casename '_det0.m']);
for i=1:length(SYS.IDX.MAT.inFlux)
    MAT(SYS.IDX.MAT.inFlux(i)).intFlux=DETintFlux(i,11);
end
SYS.intFlux=sum([MAT(SYS.IDX.MAT.inFlux).intFlux]);

%%% Get one-group XS/Nubar/Serpent k-eff/inf
run([SYS.Casename '_res.m']);
idxUni=find(all(ismember(GC_UNIVERSE_NAME,'0'),2));
idxUni=idxUni(1);
SYS.RR.NU{2}=NUBAR(idxUni,1); SYS.RR.LEAK{2} =ABS_KINF(idxUni,1)/ABS_KEFF(idxUni,1);
SYS.KEFF.Serpent(end+1)=ABS_KEFF(idxUni,1); 
SYS.KINF.Serpent(end+1)=ABS_KINF(idxUni,1);
SYS.tgtFissRate=TOT_FISSRATE(idxUni,1);

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
SYS.RR.notInMat{2}=struct('ZAI',num2cell([-1;ZAI(RROidx)]),'capt',num2cell([0;rr_sorted(RROidx,1)]),...
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
NPhi=sum(NPhi,2);
RR=RR./repmat(NPhi,1,4); % create cross-sections using integral flux and nuclides densities
RR(isnan(RR)|isinf(RR))=0.0;
SYS.RR.inMat{2}=struct('capt',RR(:,1),'fiss',RR(:,2),'n2n',RR(:,3),'n3n',RR(:,4));

%%% Read data for burnable materials
for i=SYS.IDX.MAT.burn
    run(strtrim(ls(['depmtx_' MAT(i).name '*0.m'])));
    A=sparse(A);
    idx=ismember(ZAI,[-1 10 666]);
    A(:,idx)=[]; A(idx,:)=[];
    N0(idx)=[]; ZAI(idx)=[];
    SYS.MTX.decay{2,i}=sparse(SYS.MTX.defaultDecay);
    SYS.IDX.burnZAI{2,i}=ismember(MAT(i).ZAI,ZAI);
    SYS.MTX.decay{2,i}(:,~SYS.IDX.burnZAI{2,i})=[];
    SYS.MTX.decay{2,i}(~SYS.IDX.burnZAI{2,i},:)=[];
    SYS.IDX.matZAI{2,i}=repmat(i,length(ZAI),1);
    if(~iscolumn(ZAI))
        SYS.burnZAI{2,i}=ZAI';
    else
        SYS.burnZAI{2,i}=ZAI;
    end
    SYS.MTX.burn{2,i}=sparse(A-SYS.MTX.decay{end,i});
    clearvars A N0 N1 t ZAI
end

%%% data for stream materials 
for i=SYS.IDX.REP.contMat
   SYS.MTX.decay{2,i}=SYS.MTX.decay{1,i}; 
   SYS.IDX.matZAI{2,i}=SYS.IDX.matZAI{1,i};
   SYS.IDX.burnZAI{2,i}=SYS.IDX.burnZAI{1,i};
   SYS.burnZAI{2,i}=SYS.burnZAI{1,i};
   SYS.MTX.burn{2,i}=SYS.MTX.burn{1,i};
end
return
end

