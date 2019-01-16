function [MAT,SYS] = reactivityControl(MAT,SYS,REA,IDX)
    %REACTIVITYCONTROL(MAT,SYS,REA,IDX) Controls reactivity by changing the concentrations of
    %selected nuclides in selected materials of the MAT vector according to
    %user options given in the REA vector
    global FID

    [keff,~] = computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK);
    reactDiff=[1E5*(1/REA.targetKeff-1/keff), zeros(1,REA.maxIter)]; % Reactivity in pcm
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
        if reactDiff(1)>0 %supercritical system
            nucIdx=IDX.downNucl; srcMat=IDX.targetMat; dstMat=IDX.feedMat;
            downFraction=REA.downFraction;
            if isempty(downFraction) % fractions given by source material
                downFraction=MAT(IDX.targetMat).aFrac(IDX.downNucl);
            end
            fraction=downFraction.*sum(MAT(IDX.targetMat).N(nucIdx,end));
        elseif reactDiff(1)<0 %subcritical system
            nucIdx=IDX.upNucl; dstMat=IDX.targetMat; srcMat=IDX.feedMat;
            upFraction=REA.upFraction;
            if isempty(upFraction)  % fractions given by feed material
                upFraction=MAT(IDX.feedMat).aFrac(IDX.upNucl);
            end
            fraction=upFraction.*sum(MAT(IDX.feedMat).N(nucIdx,end));
        end
        replFraction=REA.replFraction;
        if isempty(replFraction)
            replFraction=MAT(dstMat).aFrac(IDX.replNuc);
        end
        minBound=0; maxBound=1; changeUp(2)=0.1*maxBound;

        MAT(IDX.targetMat).N(:,end+1)=MAT(IDX.targetMat).N(:,end); %save copy of mat comp
        saveTarget=MAT(IDX.targetMat).N(:,end);
        if ~isempty(REA.feedMat)
            MAT(IDX.feedMat).N(:,end+1)=MAT(IDX.feedMat).N(:,end);
            saveFeed=MAT(IDX.feedMat).N(:,end);
        end

        j=1; a=minBound; b=maxBound; fa=reactDiff(1);
        while j<REA.maxIter&&abs(reactDiff(j))>REA.tol % Search loop
            j=j+1;
            MAT(IDX.targetMat).N(:,end)=saveTarget;
            if ~isempty(REA.feedMat)
                MAT(IDX.feedMat).N(:,end)=saveFeed;
            end
            if j>2
                if j==REA.maxIter
                    changeUp(j)=changeUp(max(find(abs(reactDiff)==min(abs(reactDiff)),1,'first')));
                else
                    %disp(['Search boundaries a=' num2str(a) ' and b=' num2str(b)])
                    changeUp(j)=max([a min([changeUp(j-1)-reactDiff(j-1)*(b-a)/(fb-fa) b])]); %regula falsi
                    %changeUp(j)=max([min([(changeUp(j-1)*reactDiff(j)-changeUp(j)*reactDiff(j-1))/(reactDiff(j)-reactDiff(j-1)) maxBound]),minBound]);
                end
            end
            NChange=changeUp(j)*fraction;
            switch REA.mode
            case 'addMass'
                [MAT(srcMat),MAT(dstMat)] = transferNuclides(MAT(srcMat),MAT(dstMat),nucIdx,NChange);
            case 'addVolume'
                [MAT(srcMat),MAT(dstMat)] = transferNuclides(MAT(srcMat),MAT(dstMat),nucIdx,NChange);
                MAT(srcMat).volume=MAT(srcMat).volume+changeUp(j);
            case 'replace'
                [MAT(srcMat),MAT(dstMat),NRepl] = replaceNuclides(MAT(srcMat),MAT(dstMat),nucIdx,NChange,IDX.replNuc,replFraction);
            end
            [MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat] = renormalizeSystem(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.tgtFissRate);
            [keff,~] = computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK);
            reactDiff(j)=1E5*(1/REA.targetKeff-1/keff);
            if j==2
                fb=reactDiff(2);
            else
                if reactDiff(j)*fb<0
                    b=changeUp(j); fb=reactDiff(j);
                elseif reactDiff(j)*fa<0
                    a=changeUp(j); fa=reactDiff(j);
                end
            end
            fprintf('%s\n',['  ** REACT **   Current k-eff: ' num2str(keff)]);
        end
        if reactDiff(1)>0
            NChange=NChange.*MAT(srcMat).atomicMass(IDX.downNucl)*1E24; %convert to mass
        elseif reactDiff(1)<0
            NChange=NChange.*MAT(srcMat).atomicMass(IDX.upNucl)*1E24;
        end

        if ~SYS.RUN.PCC.active||SYS.RUN.PCC.corrector % write separate output
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
                NRepl=-NRepl.*MAT(srcMat).atomicMass(IDX.replNuc)*1E24;
                NChangeWrite(ismember(elemIdx,IDX.replNuc))=NRepl;
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
