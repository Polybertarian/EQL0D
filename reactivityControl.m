function [MAT,SYS] = reactivityControl(MAT,OPT,SYS)
%REACTIVITYCONTROL(MAT,OPT,SYS) Controls reactivity by changing the concentrations of
%selected nuclides in selected materials of the MAT vector according to
%user options given in the OPT.REA vector

diff=10e5*(1/OPT.REA.targetKeff-1/SYS.keff(end)); % reactivity in pcm

if(OPT.REA.allowNegative)
    criterion=abs(diff)>OPT.REA.tol;
else
    criterion=-diff>OPT.REA.tol;
end

if(criterion) %activate reactivity control if above tolerance
    fprintf(SYS.FID.log,'%s\n',['** REACT ** k-eff: ' num2str(SYS.keff(end)) ' > tol. limit, searching...']);
    saveTarget=MAT(SYS.IDX.targetMat).N;
    if(~isempty(OPT.REA.feedMat))
        saveFeed=MAT(SYS.IDX.feedMat).N;
    end
    switch OPT.REA.mode
        case {'replace','addMass'}
            j=0;
            changeUp=[0 0];
            if(diff>0)
                maxBound=0; minBound=-1;
                changeUp(end+1)=-0.001; % initial guess
            elseif(diff<0)
                maxBound=1; minBound=0;
                changeUp(end+1)=0.001;
            end
            upEmpty=false;
            if(isempty(OPT.REA.upFraction)) % mass fractions given by feed material
                upEmpty=true;
                OPT.REA.upFraction=MAT(SYS.IDX.feedMat).mFrac(SYS.IDX.feedNucUp);
            end
            downEmpty=false;
            if(isempty(OPT.REA.downFraction))
                downEmpty=true;
                if(strcmp(OPT.REA.mode,'replace'))
                    OPT.REA.downFraction=MAT(SYS.IDX.feedMat).mFrac(SYS.IDX.feedNucDo);
                elseif(strcmp(OPT.REA.mode,'addMass'))
                    OPT.REA.downFraction=MAT(SYS.IDX.targetMat).mFrac(SYS.IDX.feedNucDo);
                end
            end
            
            while(j<OPT.REA.maxIter&abs(diff(end))>OPT.REA.tol&any(abs(changeUp(end-1:end))<1))
                j=j+1;
                if(j>1)
                    MAT(SYS.IDX.targetMat).N=saveTarget;
                    if(~isempty(OPT.REA.feedMat))
                        MAT(SYS.IDX.feedMat).N=saveFeed;
                    end
                    changeUp(end+1)=max([min([(changeUp(end-1)*diff(end)-changeUp(end)*diff(end-1))/...
                        (diff(end)-diff(end-1)) maxBound]),minBound]);
                end
                NChangeUp=0;
                NChangeDo=0;
                if(changeUp(end)<0)
                    if(strcmp(OPT.REA.mode,'replace'))
                        NChangeUp=changeUp(end)*MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucUp);
                        MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucUp)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucUp)+NChangeUp;
                        NChangeDo=OPT.REA.downFraction.*sum(NChangeUp.*MAT(SYS.IDX.targetMat).atomicMass(SYS.IDX.targetNucUp))...
                            ./MAT(SYS.IDX.targetMat).atomicMass(SYS.IDX.targetNucDo);
                    elseif(strcmp(OPT.REA.mode,'addMass'))
                        NChangeDo=changeUp(end)*MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucDo);
                        MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucDo)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucDo)+NChangeDo;
                    end
                elseif(changeUp(end)>0)
                    NChangeUp=changeUp(end)*MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucUp);
                    MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucUp)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucUp)+NChangeUp;
                    if(strcmp(OPT.REA.mode,'replace'))
                        NChangeDo=OPT.REA.downFraction.*sum(NChangeUp.*MAT(SYS.IDX.feedMat).atomicMass(SYS.IDX.feedNucUp))...
                            ./MAT(SYS.IDX.feedMat).atomicMass(SYS.IDX.feedNucDo);
                    end
                end
                if(strcmp(OPT.REA.mode,'replace'))
                    MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucDo)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucDo)-NChangeDo;
                end
                if(~isempty(OPT.REA.feedMat))
                    if(strcmp(OPT.REA.mode,'replace'))
                        MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucUp)=MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucUp)-NChangeUp;
                        MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucDo)=MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucDo)+NChangeDo;
                    elseif(strcmp(OPT.REA.mode,'addMass'))
                        if(changeUp(end)<0)
                            MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucDo)=MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucDo)-NChangeDo;
                        elseif(changeUp(end)>0)
                            MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucUp)=MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucUp)-NChangeUp;
                        end
                    end
                end
                SYS=computeK(MAT,SYS);
                diff(end+1)=10e5*(1/OPT.REA.targetKeff-1/SYS.keff(end));
                if(SYS.verboseMode)
                    fprintf(SYS.FID.log,'%s\n',['** REACT ** Current k-eff: ' num2str(SYS.keff(end))]);
                end
            end
            NChangeUp=NChangeUp.*MAT(SYS.IDX.targetMat).atomicMass(SYS.IDX.feedNucUp)*1E24;
            if(strcmp(OPT.REA.mode,'addMass'))
                if(~isempty(SYS.IDX.feedMat))
                    SYS.reactAdditions=[NChangeUp' -NChangeUp']/1000;
                else
                    SYS.reactAdditions=NChangeUp'/1000;
                end
            elseif(strcmp(OPT.REA.mode,'replace'))
                NChangeDo=NChangeDo.*MAT(SYS.IDX.targetMat).atomicMass(SYS.IDX.feedNucDo)*1E24;
                if(~isempty(SYS.IDX.feedMat))
                    SYS.reactAdditions=[NChangeUp' NChangeDo' -NChangeUp' -NChangeDo']/1000;
                else
                    SYS.reactAdditions=[NChangeUp' NChangeDo']/1000;
                end
            end
            if(upEmpty)
                OPT.REA.upFraction=[];
            end
            if(downEmpty)
                OPT.REA.downFraction=[];
            end
            if(abs(changeUp)==1)
                fprintf(SYS.FID.log,'%s\n','** REACT ** Max. composition change reached!');
            end
            if(~OPT.PCC||SYS.PCC.active)
                fprintf(SYS.FID.log,'%s\n',['** REACT ** Change: Up ' num2str(sum(NChangeUp/1000.0)) ' kg.']);
                fprintf(SYS.FID.react,'%-7.3d%-6.3d%-9G',SYS.ouCntr,SYS.inCntr,SYS.nowTime(end));
                fprintf(SYS.FID.react,[repmat('%#-13.6G',1,numel(SYS.reactAdditions)) '\n'],SYS.reactAdditions);
            end
        case 'addVolume'
            j=0;
            saveVol=MAT(SYS.IDX.targetMat).volume;
            saveTarget=MAT(SYS.IDX.targetMat).N;
            addVolume=[];
            changeVol=[0 0];
            if(diff>0)
                maxBound=0;
                minBound=-1;
                changeVol(end+1)=-0.001;
            elseif(diff<0)
                maxBound=1;
                minBound=0;
                changeVol(end+1)=0.001;
            end
            while(j<OPT.REA.maxIter&abs(diff(end))>OPT.REA.tol)
                j=j+1;
                if(j>1)
                    MAT(SYS.IDX.targetMat).N=saveTarget;
                    MAT(SYS.IDX.targetMat).volume=saveVol;
                    changeVol(end+1)=max([min([(changeVol(end-1)*diff(end)-changeVol(end)*diff(end-1))/...
                        (diff(end)-diff(end-1)) maxBound]), minBound]);
                end
                addVolume(end+1)=changeVol(end)*saveVol;
                MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNuc)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNuc)+...
                    addVolume(end)*MAT(SYS.IDX.feedMat).atDens(SYS.IDX.feedNuc);
                MAT(SYS.IDX.targetMat).volume=MAT(SYS.IDX.targetMat).volume+addVolume(end);
                SYS=computeK(MAT,SYS);
                diff(end+1)=10e5*(1/OPT.REA.targetKeff-1/SYS.keff(end));
            end
            MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNuc)=MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNuc)-addVolume(end)*MAT(SYS.IDX.feedMat).atDens(SYS.IDX.feedNuc);
            if(SYS.verboseMode)
                fprintf(SYS.FID.log,'%s\n',['** REACT ** Current k-eff: ' num2str(SYS.keff(end))]);
            end
            if(~OPT.PCC||SYS.PCC.active)
                MAT(SYS.IDX.feedMat).volume=MAT(SYS.IDX.feedMat).volume-addVolume(end);
                fprintf(SYS.FID.volume,'%-7.3d %-6.3d %-9G %#-13.6G\n',...
                    SYS.ouCntr,SYS.inCntr,SYS.nowTime(end),addVolume(end));
            end
    end
    if(j>=OPT.REA.maxIter)
        fprintf(SYS.FID.log,'%s\n','** REACT ** Max. # iterations reached!');
    else
        fprintf(SYS.FID.log,'%s\n',['** REACT ** k-eff: ' num2str(SYS.keff(end)) ', stopping... ']);
    end
else
    fprintf(SYS.FID.log,'%s\n',['** REACT ** k-eff: ' num2str(SYS.keff(end)) ' < tol. limit, no search.']);
end
return
end

