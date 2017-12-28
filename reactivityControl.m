function [MAT,SYS] = reactivityControl(MAT,OPT,SYS)
%REACTIVITYCONTROL(MAT,OPT,SYS) Controls reactivity by changing the concentrations of
%selected nuclides in selected materials of the MAT vector according to
%user options given in the OPT.REA vector

diff=1e5*(1/OPT.REA.targetKeff-1/SYS.keff(end)); % reactivity in pcm

if(OPT.REA.allowNegative)
    criterion=abs(diff)>OPT.REA.tol;
else
    criterion=diff<-OPT.REA.tol;
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
            changeUp=0;
            if(diff>0)
                maxBound=0; minBound=-1;
                changeUp(end+1)=-0.99; % initial guess
                replReset=false;
                if(isempty(OPT.REA.replFraction)&strcmp(OPT.REA.mode,'replace'))
                  replReset=true;
                  OPT.REA.replFraction=MAT(SYS.IDX.feedMat).mFrac(SYS.IDX.feedNucRepl);
                end
            elseif(diff<0)
                maxBound=1; minBound=0;
                changeUp(end+1)=0.99;
                replReset=false;
                if(isempty(OPT.REA.replFraction)&strcmp(OPT.REA.mode,'replace'))
                  replReset=true;
                  OPT.REA.replFraction=MAT(SYS.IDX.targetMat).mFrac(SYS.IDX.targetNucRepl);
                end
            end
            upReset=false;
            if(isempty(OPT.REA.upFraction)) % mass fractions given by feed material
                upReset=true;
                OPT.REA.upFraction=MAT(SYS.IDX.feedMat).mFrac(SYS.IDX.feedNucUp);
            end
            downReset=false;
            if(isempty(OPT.REA.downFraction))
                downReset=true;
                OPT.REA.downFraction=MAT(SYS.IDX.targetMat).mFrac(SYS.IDX.targetNucDo);
            end
            while(j<OPT.REA.maxIter&abs(diff(end))>OPT.REA.tol)
                j=j+1;
                if(j>1)
                    MAT(SYS.IDX.targetMat).N=saveTarget;
                    if(~isempty(OPT.REA.feedMat))
                        MAT(SYS.IDX.feedMat).N=saveFeed;
                    end
                    if(j==OPT.REA.maxIter)
                        changeUp(end+1)=changeUp(abs(diff)==min(abs(diff)));
                    else
                    	changeUp(end+1)=max([min([(changeUp(end-1)*diff(end)-changeUp(end)*diff(end-1))/(diff(end)-diff(end-1)) maxBound]),minBound]);
                    end
                end
                NChange=0;
                if(diff(1)>0)
                  NChange=changeUp(end)*MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucDo);
                  MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucDo)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucDo)+NChange;
                  if(strcmp(OPT.REA.mode,'replace'))
                    NChangeRepl=OPT.REA.replFraction.*sum(NChange.*MAT(SYS.IDX.targetMat).atomicMass(SYS.IDX.targetNucDo))...
                        ./MAT(SYS.IDX.targetMat).avMass(SYS.IDX.targetNucRepl);
                    MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucRepl)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucRepl)-NChangeRepl;
                  end
                elseif(diff(1)<0)
                  NChange=changeUp(end)*MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucUp);
                  MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucUp)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucUp)+NChange;
                  if(strcmp(OPT.REA.mode,'replace'))
                    NChangeRepl=OPT.REA.replFraction.*sum(NChange.*MAT(SYS.IDX.feedMat).atomicMass(SYS.IDX.feedNucUp))...
                        ./MAT(SYS.IDX.feedMat).avMass(SYS.IDX.feedNucRepl);
                    MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucRepl)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNucRepl)-NChangeRepl;
                  end
                end
                if(~isempty(OPT.REA.feedMat))
                  if(diff(1)>0)
                      MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucDo)=MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucDo)-NChange;
                  elseif(diff(1)<0)
                      MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucUp)=MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucUp)-NChange;
                  end
                  if(strcmp(OPT.REA.mode,'replace'))
                    MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucRepl)=MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNucRepl)+NChangeRepl;
                  end
                end
                SYS=computeK(MAT,SYS);
                diff(end+1)=1e5*(1/SYS.keff(end)-1/OPT.REA.targetKeff);
                if(SYS.verboseMode)
                    fprintf(SYS.FID.log,'%s\n',['** REACT ** Current k-eff: ' num2str(SYS.keff(end))]);
                end
            end
            if(diff>0)
              NChange=NChange.*MAT(SYS.IDX.targetMat).atomicMass(SYS.IDX.targetNucDo)*1E24;
            elseif(diff<0)
              NChange=NChange.*MAT(SYS.IDX.targetMat).atomicMass(SYS.IDX.targetNucUp)*1E24;
            end
            if(strcmp(OPT.REA.mode,'addMass'))
                if(~isempty(SYS.IDX.feedMat))
                    SYS.reactAdditions=[NChange' -NChange']/1000;
                else
                    SYS.reactAdditions=NChange'/1000;
                end
            elseif(strcmp(OPT.REA.mode,'replace'))
                NChangeRepl=NChangeRepl.*MAT(SYS.IDX.targetMat).atomicMass(SYS.IDX.feedNucRepl)*1E24;
                if(~isempty(SYS.IDX.feedMat))
                    SYS.reactAdditions=[NChange' -NChangeRepl' -NChange' +NChangeRepl']/1000;
                else
                    SYS.reactAdditions=[NChange' -NChangeRepl']/1000;
                end
            end
            if(upReset)
                OPT.REA.upFraction=[];
            end
            if(downReset)
                OPT.REA.downFraction=[];
            end
            if(replReset)
                OPT.REA.replFraction=[];
            end
            if(abs(changeUp)==1)
                fprintf(SYS.FID.log,'%s\n','** REACT ** Max. composition change reached!');
            end
            if(~OPT.PCC||SYS.PCC.corrector)
                fprintf(SYS.FID.log,'%s\n',['** REACT ** Change: ' num2str(sum(NChange/1000.0)) ' kg.']);
                fprintf(SYS.FID.react,'%-7.3d%-6.3d%-9G',SYS.ouCntr,SYS.inCntr,SYS.nowTime(end));
                fprintf(SYS.FID.react,[repmat('%#-13.6G',1,numel(SYS.reactAdditions)) '\n'],SYS.reactAdditions);
            end
        case 'addVolume'
            j=0;
            saveVol=MAT(SYS.IDX.targetMat).volume;
            saveTarget=MAT(SYS.IDX.targetMat).N;
            addVolume=[];
            changeVol=0;
            if(diff<0)
                maxBound=0;
                minBound=-1;
                changeVol(end+1)=-0.05;
            elseif(diff>0)
                maxBound=1;
                minBound=0;
                changeVol(end+1)=0.05;
            end
            while(j<OPT.REA.maxIter&abs(diff(end))>OPT.REA.tol)
                j=j+1;
                if(j>1)
                    MAT(SYS.IDX.targetMat).N=saveTarget;
                    MAT(SYS.IDX.targetMat).volume=saveVol;
                    if(j==OPT.REA.maxIter)
                        changeVol(end+1)=changeVol(abs(diff)==min(abs(diff)));
                    else
                        changeVol(end+1)=max([min([(changeVol(end-1)*diff(end)-changeVol(end)*diff(end-1))/(diff(end)-diff(end-1)) maxBound]), minBound]);
                    end
                end
                addVolume(end+1)=changeVol(end)*saveVol;
                MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNuc)=MAT(SYS.IDX.targetMat).N(SYS.IDX.targetNuc)+addVolume(end)*MAT(SYS.IDX.feedMat).atDens(SYS.IDX.feedNuc);
                MAT(SYS.IDX.targetMat).volume=MAT(SYS.IDX.targetMat).volume+addVolume(end);
                SYS=computeK(MAT,SYS);
                diff(end+1)=1e5*(1/SYS.keff(end)-1/OPT.REA.targetKeff);
            end
            MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNuc)=MAT(SYS.IDX.feedMat).N(SYS.IDX.feedNuc)-addVolume(end)*MAT(SYS.IDX.feedMat).atDens(SYS.IDX.feedNuc);
            if(SYS.verboseMode)
                fprintf(SYS.FID.log,'%s\n',['** REACT ** Current k-eff: ' num2str(SYS.keff(end))]);
            end
            if(~OPT.PCC||SYS.PCC.corrector)
                MAT(SYS.IDX.feedMat).volume=MAT(SYS.IDX.feedMat).volume-addVolume(end);
                fprintf(SYS.FID.volume,'%-7.3d %-6.3d %-9G %#-13.6G\n',SYS.ouCntr,SYS.inCntr,SYS.nowTime(end),addVolume(end));
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
