function [MAT,SYS] = burnCycle(MAT,OPT,REP,SYS)
%[MAT,SYS] = burnCycle(MAT,OPT,REP,SYS) depletes the materials in the
%System
%% Burning
if(SYS.PCC.corrector)
    prefix='C';
    prefix2='CORRECTOR';
    if(SYS.inCntr==1)
        SYS.nowTime(end+1)=SYS.nowTime(end)-SYS.tStep(end)*OPT.nSteps(SYS.ouCntr)/(24.0*3600.0);
        for i=SYS.IDX.contMat
            MAT(i).N(:,end+1)=SYS.prevN.BOC(:,i);
        end
    end
    MAT=updateRates(MAT,SYS);
    SYS.MTX.interp=((SYS.PCC.nSteps-SYS.inCntr+1)*SYS.MTX.total{1}+...
        (SYS.inCntr-1)*SYS.MTX.total{2})/SYS.PCC.nSteps; %%% Interpolate and solve
    SYS.burnVec=CRAMsolve(SYS.MTX.interp,SYS.tStep(end),SYS.burnVec);
    SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0);
    for i=SYS.IDX.contMat %%% Change material compositions
        MAT(i).N(SYS.IDX.burnZAI{2,i},end+1)=SYS.burnVec(vertcat(SYS.IDX.matZAI{2,:})==i)*MAT(i).volume;
    end
else
    prefix='P';
    prefix2='PREDICTOR';
    SYS.burnVec=[]; %%% Prepare vector for solving
    for i=SYS.IDX.contMat
        SYS.burnVec=vertcat(SYS.burnVec,MAT(i).atDens(SYS.IDX.burnZAI{2,i}));
    end
    if(SYS.PCC.active)
        SYS.burnVec2=CRAMsolve(SYS.MTX.total{2},SYS.tStep(end)*OPT.nSteps(SYS.ouCntr),SYS.burnVec);
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)*OPT.nSteps(SYS.ouCntr)/(24.0*3600.0); %in EFPD
    else
        SYS.burnVec2=CRAMsolve(SYS.MTX.total{2},SYS.tStep(end),SYS.burnVec);
        SYS.nowTime(end+1)=SYS.nowTime(end)+SYS.tStep(end)/(24.0*3600.0); %in EFPD
    end
    for i=SYS.IDX.contMat %%% Change material compositions
        MAT(i).N(SYS.IDX.burnZAI{2,i},end+1)=SYS.burnVec2(vertcat(SYS.IDX.matZAI{2,:})==i)*MAT(i).volume;
    end
end

%% Before Batch Processing

if(strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10))
    printStatus(MAT,SYS,prefix2,'EoS (BB)');            %%% Add status to log
end

SYS=computeK(MAT,SYS); %%% Compute k-eff
printK(SYS,'BB',prefix,'EQL0D');

if(OPT.printSteps&&OPT.printStepsBB)
    for i=SYS.IDX.contMat
        MAT(i).printMaterial(SYS,'BB'); %%% Print EOS composition
    end
end

SYS.prevN.EOC=[];
for i=1:length(MAT)
    SYS.prevN.EOC=[SYS.prevN.EOC MAT(i).N(:,end)]; %%% Store compositions from previous loop
end

if(OPT.redoxControl) %%% Adjust redox
    MAT = redoxControl(MAT,OPT,SYS);
end
if(~isempty(SYS.IDX.batchStr)) %%% Batchwise EoS processes
    MAT=batchProcessing(MAT,REP,SYS);
end
SYS=computeK(MAT,SYS); %%% Compute k-eff
if(OPT.reactControl) %%% Adjust reactivity
    [MAT,SYS] = reactivityControl(MAT,OPT,SYS);
end

%% After Batch processing
if(strcmp(OPT.iterMode,'steps')||SYS.inCntr==1||floor(SYS.inCntr/10)==ceil(SYS.inCntr/10))
    printStatus(MAT,SYS,prefix2,'EoS (AB)');  %%% Add status to log
end

if(SYS.inCntr<OPT.nSteps(SYS.ouCntr))
    SYS=computeK(MAT,SYS); %%% Compute k-eff
    printK(SYS,'AB',prefix,'EQL0D');
end

%%% Print material composition to file
if(OPT.printSteps)
    for i=SYS.IDX.contMat
        MAT(i).printMaterial(SYS,'AB');
    end
end

for i=SYS.IDX.burnMat
    MAT(i).write(OPT.matWriteStyle); %%% Write compositions
end

end

