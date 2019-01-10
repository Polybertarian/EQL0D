function [MAT,SYS] = reactivityControl(MAT,SYS,REA,IDX)
    %REACTIVITYCONTROL(MAT,SYS,REA,IDX) Controls reactivity by changing the concentrations of
    %selected nuclides in selected materials of the MAT vector according to
    %user options given in the REA vector
    global FID

    if(any(isnan(MAT(IDX.feedMat).N(:,end))))
        error('lol')
    end

    reactDiff=[1E5*(1/REA.targetKeff-1/SYS.KEFF.EQL0D(end)), zeros(1,REA.maxIter)]; % Reactivity in pcm
    if REA.allowRemoval&&REA.allowAddition
        criterion=abs(reactDiff(1))>abs(REA.tol);
    elseif REA.allowRemoval&&~REA.allowAddition
        criterion=reactDiff(1)>abs(REA.tol);
    elseif REA.allowAddition&&~REA.allowRemoval
        criterion=reactDiff(1)<-abs(REA.tol);
    end

    if ~criterion % Activate reactivity control if above tolerance
        fprintf('%s\n',['  ** REACT **   k-eff: ' num2str(SYS.KEFF.EQL0D(end)) ' < tol. limit, no search.']);
    else
        fprintf('%s\n',['  ** REACT **   k-eff: ' num2str(SYS.KEFF.EQL0D(end)) ' > tol. limit, searching...']);
        changeUp=zeros(1,REA.maxIter+1);
        upFraction=REA.upFraction; downFraction=REA.downFraction; replFraction=REA.replFraction;
        if isempty(upFraction)  % Mass fractions given by feed material
            upFraction=MAT(IDX.feedMat).mFrac(IDX.upNucl);
        end
        if isempty(downFraction) % Mass fractions given by source material
            downFraction=MAT(IDX.targetMat).mFrac(IDX.downNucl);
        end
        if reactDiff(1)>0
            maxBound=0; minBound=-1; changeUp(2)=-0.1;
            if isempty(replFraction)
                replFraction=MAT(IDX.feedMat).mFrac(IDX.replNuc);
            end
        elseif reactDiff(1)<0
            maxBound=1; minBound=0; changeUp(2)=0.1;
            if isempty(replFraction)
                replFraction=MAT(IDX.targetMat).mFrac(IDX.replNuc);
            end
        end

        %% Search loop
        j=0;
        while j<REA.maxIter&&abs(reactDiff(j+1))>REA.tol
            j=j+1;
            if j==1
                MAT(IDX.targetMat).N(:,end+1)=MAT(IDX.targetMat).N(:,end);
                if ~isempty(REA.feedMat)
                    MAT(IDX.feedMat).N(:,end+1)=MAT(IDX.feedMat).N(:,end);
                end
            elseif j>1
                MAT(IDX.targetMat).N(:,end)=MAT(IDX.targetMat).N(:,end-1);
                if ~isempty(REA.feedMat)
                    MAT(IDX.feedMat).N(:,end)=MAT(IDX.feedMat).N(:,end-1);
                end
                if j==REA.maxIter
                    changeUp(j+1)=changeUp(max(find(abs(reactDiff)==min(abs(reactDiff)),1,'first')));
                else
                    changeUp(j+1)=max([min([(changeUp(j-1)*reactDiff(j)-changeUp(j)*reactDiff(j-1))/...
                    (reactDiff(j)-reactDiff(j-1)) maxBound]),minBound]);
                end
            end
            NChange=0;
            if reactDiff(1)>0
                NChange=changeUp(j+1)*downFraction.*sum(MAT(IDX.targetMat).N(IDX.downNucl,end));
                switch REA.mode
                case 'addMass'
                    [MAT(IDX.targetMat),MAT(IDX.feedMat)] = transferNuclides(MAT(IDX.targetMat),MAT(IDX.feedMat),IDX.downNucl,NChange);
                case 'addVolume'
                    [MAT(IDX.targetMat),MAT(IDX.feedMat)] = transferNuclides(MAT(IDX.targetMat),MAT(IDX.feedMat),IDX.upNucl,NChange);
                    MAT(IDX.targetMat).volume=MAT(IDX.targetMat).volume+changeUp(j);
                case 'replace'
                    [MAT(IDX.targetMat),MAT(IDX.feedMat),NChangeRepl] = replaceNuclides(MAT(IDX.targetMat),MAT(IDX.feedMat),...
                    IDX.downNucl,NChange,IDX.replNuc,replFraction);
                end
            elseif reactDiff(1)<0
                NChange=changeUp(j+1)*upFraction.*sum(MAT(IDX.feedMat).N(IDX.upNucl,end));
                switch REA.mode
                case 'addMass'
                    [MAT(IDX.targetMat),MAT(IDX.feedMat)] = transferNuclides(MAT(IDX.targetMat),MAT(IDX.feedMat),IDX.upNucl,NChange);
                case 'addVolume'
                    [MAT(IDX.targetMat),MAT(IDX.feedMat)] = transferNuclides(MAT(IDX.targetMat),MAT(IDX.feedMat),IDX.upNucl,NChange);
                    MAT(IDX.targetMat).volume=MAT(IDX.targetMat).volume+changeUp(j);
                case 'replace'
                    [MAT(IDX.targetMat),MAT(IDX.feedMat),NChangeRepl] = replaceNuclides(MAT(IDX.targetMat),MAT(IDX.feedMat),...
                    IDX.upNucl,NChange,IDX.replNuc,replFraction);
                end
            end
            [MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat] = renormalizeSystem(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.tgtFissRate);
            [keff,~] = computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK);
            reactDiff(j+1)=1e5*(1/keff-1/REA.targetKeff);
            fprintf('%s\n',['  ** REACT **   Current k-eff: ' num2str(SYS.KEFF.EQL0D(end))]);
        end
        if reactDiff(1)>0
            NChange=NChange.*MAT(IDX.targetMat).atomicMass(IDX.downNucl)*1E24;
        elseif reactDiff(1)<0
            NChange=NChange.*MAT(IDX.targetMat).atomicMass(IDX.upNucl)*1E24;
        end
        if ~SYS.RUN.PCC.active||SYS.RUN.PCC.corrector
            if strcmp(REA.mode,'addMass')
                elemIdx=unique([IDX.upNucl;IDX.downNucl]);
                NChangeWrite=zeros(size(elemIdx));
                if reactDiff(1)>0
                    NChangeWrite(ismember(elemIdx,IDX.downNucl))=NChange;
                elseif reactDiff(1)<0
                    NChangeWrite(ismember(elemIdx,IDX.upNucl))=NChange;
                end
            elseif strcmp(REA.mode,'replace')
                elemIdx=unique([IDX.replNuc;IDX.upNucl;IDX.downNucl]);
                NChangeWrite=zeros(size(elemIdx));
                if reactDiff(1)>0
                    NChangeWrite(ismember(elemIdx,IDX.downNucl))=NChange;
                elseif reactDiff(1)<0
                    NChangeWrite(ismember(elemIdx,IDX.upNucl))=NChange;
                end
                NChangeRepl=-NChangeRepl.*MAT(IDX.targetMat).atomicMass(IDX.replNuc)*1E24;
                NChangeWrite(ismember(elemIdx,IDX.replNuc))=NChangeRepl;
            end
            if ~isempty(IDX.feedMat)
                REActAdditions=[NChangeWrite' -NChangeWrite']/1000;
            else
                REActAdditions=NChangeWrite'/1000;
            end
            if abs(changeUp(j))==maxBound||abs(changeUp(j))==minBound
                fprintf('%s\n','  ** REACT **   Max. composition change reached!');
            end
            switch REA.mode
            case {'addMass','replace'}
                fprintf('%s\n',['  ** REACT **   Time: ' num2str(SYS.RUN.nowTime(end)) ...
                ' EFPD Change: ' num2str(sum(NChange/1000.0),'%8.6G') ' kg.']);
                fprintf(FID.react,'%-7.3d%-6.3d%-9G',SYS.RUN.ouCntr,SYS.RUN.inCntr,SYS.RUN.nowTime(end));
                fprintf(FID.react,[repmat('%-13.6G',1,numel(REActAdditions)) '\n'],REActAdditions);
            case 'addVolume'
                MAT(IDX.feedMat).volume=MAT(IDX.feedMat).volume-addVolume(j);
                fprintf(FID.volume,'%-7.3d %-6.3d %-9G %#-13.6G\n',SYS.RUN.ouCntr,SYS.RUN.inCntr,SYS.RUN.nowTime(end),addVolume(j));
            end
            if j>=REA.maxIter
                fprintf('%s\n','  ** REACT **   Max. # iterations reached!');
            else
                fprintf('%s\n',['  ** REACT **   k-eff: ' num2str(SYS.KEFF.EQL0D(end)) ', stopping... ']);
            end
        end
    end
