function [MAT,SYS] = reactivityControl(MAT,SYS,REA,IDX)
    %REACTIVITYCONTROL(MAT,SYS,REA,IDX) Controls reactivity by changing the concentrations of
    %selected nuclides in selected materials of the MAT vector according to
    %user options given in the REA vector

    global FID

    reactDiff=[1E5*(1/REA.targetKeff-1/SYS.KEFF.EQL0D(end)), zeros(1,REA.maxIter)]; % Reactivity in pcm

    if REA.allowRemoval&&REA.allowAddition
        criterion=abs(reactDiff(1))>abs(REA.tol);
    elseif REA.allowRemoval&&~REA.allowAddition
        criterion=reactDiff(1)>abs(REA.tol);
    elseif REA.allowAddition&&~REA.allowRemoval
        criterion=reactDiff(1)<-abs(REA.tol);
    end

    if criterion % Activate reactivity control if above tolerance
        saveTarget=MAT(IDX.target).N(:,end);
        if ~isempty(REA.feedMat)
            saveFeed=MAT(IDX.feed).N(:,end);
        end
        fprintf('%s\n',['** REACT ** k-eff: ' num2str(SYS.KEFF.EQL0D(end)) ' > tol. limit, searching...']);

        switch REA.mode
        case {'replace','addMass'}
            j=0;
            changeUp=zeros(1,REA.maxIter+1);
            if reactDiff(1)>0
                maxBound=0; minBound=-1; changeUp(j+1)=-0.99;replReset=false;
                if isempty(REA.replFraction)&&strcmp(REA.mode,'replace')
                    replReset=true;
                    REA.replFraction=MAT(IDX.feed).mFrac(IDX.feedNucRepl);
                end
            elseif reactDiff(1)<0
                maxBound=1; minBound=0; changeUp(j+1)=0.99; replReset=false;
                if isempty(REA.replFraction)&&strcmp(REA.mode,'replace')
                    replReset=true;
                    REA.replFraction=MAT(IDX.target).mFrac(IDX.targetNucRepl);
                end
            end
            upReset=false;
            if isempty(REA.upFraction)  % Mass fractions given by feed material
                upReset=true;
                REA.upFraction=MAT(IDX.feed).mFrac(IDX.feedNucUp);
            end
            downReset=false;
            if isempty(REA.downFraction)
                downReset=true;
                REA.downFraction=MAT(IDX.target).mFrac(IDX.targetNucDo);
            end
            while j<REA.maxIter&&abs(reactDiff(j+1))>REA.tol
                j=j+1;
                if j==1
                    MAT(IDX.target).N(:,end+1)=saveTarget;
                    if ~isempty(REA.feedMat)
                        MAT(IDX.feed).N(:,end+1)=saveFeed;
                    end
                elseif j>1
                    MAT(IDX.target).N(:,end)=saveTarget;
                    if ~isempty(REA.feedMat)
                        MAT(IDX.feed).N(:,end)=saveFeed;
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
                    NChange=changeUp(j)*MAT(IDX.target).N(IDX.targetNucDo,end);
                    MAT(IDX.target).N(IDX.targetNucDo,end)=MAT(IDX.target).N(IDX.targetNucDo,end)+NChange;
                    if strcmp(REA.mode,'replace')
                        NChangeRepl=REA.replFraction.*sum(NChange.*MAT(IDX.target).atomicMass(IDX.targetNucDo))...
                        ./MAT(IDX.target).avMass(IDX.targetNucRepl);
                        MAT(IDX.target).N(IDX.targetNucRepl,end)=MAT(IDX.target).N(IDX.targetNucRepl,end)-NChangeRepl;
                    end
                elseif reactDiff(1)<0
                    NChange=changeUp(j)*MAT(IDX.feed).N(IDX.feedNucUp,end);
                    MAT(IDX.target).N(IDX.targetNucUp,end)=MAT(IDX.target).N(IDX.targetNucUp,end)+NChange;
                    if strcmp(REA.mode,'replace')
                        NChangeRepl=REA.replFraction.*sum(NChange.*MAT(IDX.feed).atomicMass(IDX.feedNucUp))...
                        ./MAT(IDX.feed).avMass(IDX.feedNucRepl);
                        MAT(IDX.target).N(IDX.targetNucRepl,end)=MAT(IDX.target).N(IDX.targetNucRepl,end)-NChangeRepl;
                    end
                end
                if ~isempty(REA.feedMat)
                    if reactDiff(1)>0
                        MAT(IDX.feed).N(IDX.feedNucDo,end)=MAT(IDX.feed).N(IDX.feedNucDo,end)-NChange;
                    elseif reactDiff(1)<0
                        MAT(IDX.feed).N(IDX.feedNucUp,end)=MAT(IDX.feed).N(IDX.feedNucUp,end)-NChange;
                    end
                    if strcmp(REA.mode,'replace')
                        MAT(IDX.feed).N(IDX.feedNucRepl,end)=MAT(IDX.feed).N(IDX.feedNucRepl,end)+NChangeRepl;
                    end
                end
                [MAT(SYS.IDX.MAT.inFlux),SYS.RR] = renormalizeSystem(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.tgtFissRate)
                [keff,~]=computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK);
                reactDiff(j+1)=1e5*(1/keff-1/REA.targetKeff);
                fprintf('%s\n',['** REACT ** Current k-eff: ' num2str(SYS.KEFF.EQL0D(end))]);
            end
            if reactDiff(1)>0
                NChange=NChange.*MAT(IDX.target).atomicMass(IDX.targetNucDo)*1E24;
            elseif reactDiff(1)<0
                NChange=NChange.*MAT(IDX.target).atomicMass(IDX.targetNucUp)*1E24;
            end
            if strcmp(REA.mode,'addMass')
                elemIdx=unique([IDX.targetNucUp;IDX.targetNucDo]);
                NChangeWrite=zeros(size(elemIdx));
                if reactDiff(1)>0
                    NChangeWrite(ismember(elemIdx,IDX.targetNucDo))=NChange;
                elseif reactDiff(1)<0
                    NChangeWrite(ismember(elemIdx,IDX.targetNucUp))=NChange;
                end
            elseif strcmp(REA.mode,'replace')
                elemIdx=unique([IDX.targetNucRepl;IDX.targetNucUp;IDX.targetNucDo]);
                NChangeWrite=zeros(size(elemIdx));
                if reactDiff(1)>0
                    NChangeWrite(ismember(elemIdx,IDX.targetNucDo))=NChange;
                elseif reactDiff(1)<0
                    NChangeWrite(ismember(elemIdx,IDX.targetNucUp))=NChange;
                end
                NChangeRepl=-NChangeRepl.*MAT(IDX.target).atomicMass(IDX.feedNucRepl)*1E24;
                NChangeWrite(ismember(elemIdx,IDX.targetNucRepl))=NChangeRepl;
            end
            if ~isempty(IDX.feed)
                REActAdditions=[NChangeWrite' -NChangeWrite']/1000;
            else
                REActAdditions=NChangeWrite'/1000;
            end
            if upReset
                REA.upFraction=[];
            end
            if downReset
                REA.downFraction=[];
            end
            if replReset
                REA.replFraction=[];
            end
            if abs(changeUp)==1
                fprintf('%s\n','** REACT ** Max. composition change reached!');
            end
            if ~SYS.RUN.PCC.active||SYS.RUN.PCC.corrector
                fprintf('%s\n',['** REACT ** Time: ' num2str(SYS.RUN.nowTime(end)) ...
                ' EFPD Change: ' num2str(sum(NChange/1000.0),'%8.6G') ' kg.']);
                fprintf(FID.react,'%-7.3d%-6.3d%-9G',SYS.RUN.ouCntr,SYS.RUN.inCntr,SYS.RUN.nowTime(end));
                fprintf(FID.react,[repmat('%-13.6G',1,numel(REActAdditions)) '\n'],REActAdditions);
            end
        case 'addVolume'
            j=0;
            saveVol=MAT(IDX.target).volume;
            addVolume=zeros(1,REA.maxIter);
            changeVol=zeros(1,REA.maxIter+1);
            if reactDiff(1)<0
                maxBound=0; minBound=-1; changeVol(j+1)=-0.05;
            elseif reactDiff(1)>0
                maxBound=1; minBound=0; changeVol(j+1)=0.05;
            end
            while j<REA.maxIter&&abs(reactDiff(j))>REA.tol
                j=j+1;
                if j==1
                    MAT(IDX.target).N(:,end+1)=saveTarget;
                    if ~isempty(REA.feedMat)
                        MAT(IDX.feed).N(:,end+1)=saveFeed;
                    end
                elseif j>1
                    MAT(IDX.target).N(:,end)=saveTarget;
                    MAT(IDX.target).volume=saveVol;
                    if j==REA.maxIter
                        changeVol(j+1)=changeVol(max(find(abs(reactDiff)==min(abs(reactDiff)),1,'first')));
                    else
                        changeVol(j+1)=max([min([(changeVol(j-1)*reactDiff(j)-changeVol(j)*reactDiff(j-1))/...
                        (reactDiff(j)-reactDiff(j-1)) maxBound]), minBound]);
                    end
                end
                addVolume(j+1)=changeVol(j)*saveVol;
                MAT(IDX.target).N(IDX.targetNucUp,end)=...
                MAT(IDX.target).N(IDX.targetNucUp,end)+addVolume(j)*MAT(IDX.feed).atDens(IDX.feedNuc);
                MAT(IDX.target).volume=MAT(IDX.target).volume+addVolume(j);
                [keff,~]=computeK(MAT(SYS.IDX.MAT.inFlux),SYS.RR(3).notInMat,SYS.NUBAR,SYS.LEAK);
                reactDiff(j+1)=1e5*(1/keff-1/REA.targetKeff);
            end
            MAT(IDX.feed).N(IDX.feedNuc,end)=...
            MAT(IDX.feed).N(IDX.feedNuc,end)-addVolume(j)*MAT(IDX.feed).atDens(IDX.feedNuc);

            fprintf('%s\n',['** REACT ** Current k-eff: ' num2str(SYS.KEFF.EQL0D(end))]);
            if ~SYS.RUN.PCC.active||SYS.RUN.PCC.corrector
                MAT(IDX.feed).volume=MAT(IDX.feed).volume-addVolume(j);
                fprintf(FID.volume,'%-7.3d %-6.3d %-9G %#-13.6G\n',SYS.RUN.ouCntr,SYS.RUN.inCntr,SYS.RUN.nowTime(end),addVolume(j));
            end
        end
        if j>=REA.maxIter
            fprintf('%s\n','** REACT ** Max. # iterations reached!');
        else
            fprintf('%s\n',['** REACT ** k-eff: ' num2str(SYS.KEFF.EQL0D(end)) ', stopping... ']);
        end
    else
        fprintf('%s\n',['** REACT ** k-eff: ' num2str(SYS.KEFF.EQL0D(end)) ' < tol. limit, no search.']);
    end
    return
end
