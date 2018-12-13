function SYS = testConvergence(MAT,OPT,SYS,loop)
%TESTCONVERGENCE tests convergence in loop of Mat composition based on convParam

global FID

switch OPT.iterMode
  case 'steps'
    if(strcmp(loop,'outer'))
      SYS.stopOuter=(SYS.ouCntr>=OPT.nCycles);
      return
    else
      SYS.stopInner=(SYS.inCntr>=OPT.nSteps(SYS.ouCntr));
      return
    end
  case 'equilibrium'
    if(strcmp(loop,'outer'))
      if(SYS.ouCntr>=OPT.nCycles)
        SYS.stopOuter=true;
        fprintf('%s\n','*** OUTER CONVERGENCE *** maximum iterations reached!');
        return
      elseif(SYS.ouCntr==1)
        SYS.stopOuter=false;
        return
      end
    else
      if(SYS.inCntr>=OPT.nSteps(SYS.ouCntr))
        SYS.stopInner=true;
        fprintf('%s\n','** INNER CONVERGENCE ** maximum iterations reached!');
        return
      end
    end
    isConverged=zeros(size(SYS.IDX.MAT.burn));
    if(ismember(OPT.CONV.(loop).criterion,{'maxRelDiff','maxActRelDiff','maxNucDatRelDiff'}))
      %%% Loop over burnt materials
      for i=SYS.IDX.MAT.burn
        switch OPT.CONV.(loop).criterion
          case 'maxRelDiff'
            idx=find(MAT(i).atDens>OPT.CONV.(loop).cutoff);
          case 'maxNucDatRelDiff'
            idx=find(MAT(i).hasNucData&MAT(i).atDens>OPT.CONV.(loop).cutoff);
          case 'maxActRelDiff'
            idx=find(isActinide(MAT(i).ZAI)&MAT(i).atDens>OPT.CONV.(loop).cutoff);
        end
        
        %%% Find reference composition
        switch loop
          case 'inner'
            Nold=SYS.prevN.BOC(:,i);
          case 'outer'
            Nold=SYS.oldN(:,i);
        end
        
        %%% Compute relative difference
        relDiff=abs((Nold(idx)-MAT(i).N(idx,end))./Nold(idx));
        relDiff(isnan(relDiff))=0;
        relDiff(isinf(relDiff))=0;
        
        [diff,idx2]=max(relDiff);
        switch loop
          case 'inner'
            if(floor(SYS.inCntr/20)==ceil(SYS.inCntr/20))
              fprintf('%s\n',['** INNER CONVERGENCE ** outer ' num2str(SYS.ouCntr) ' inner ' ...
                num2str(SYS.inCntr) ': ' OPT.CONV.(loop).criterion ' in ' MAT(i).name ' for ' ...
                char(ZAI2Name(MAT(i).ZAI(idx(idx2)))) ': ' num2str(diff,'%.3G')]);
            end
          case 'outer'
            fprintf('%s\n',['*** OUTER CONVERGENCE *** outer ' num2str(SYS.ouCntr) ': ' ...
              OPT.CONV.(loop).criterion ' in ' MAT(i).name ' for ' char(ZAI2Name(MAT(i).ZAI(idx(idx2)))) ': ' ...
              num2str(diff,'%.3G')]);
        end
        if(diff<OPT.CONV.(loop).value)
          isConverged(i)=true;
        else
          isConverged(i)=false;
        end
        
      end
    else
      for i=SYS.IDX.MAT.burn
        switch OPT.CONV.(loop).criterion
          case 'maxBUDiff'
            switch loop
              case 'inner'
                diff=(MAT(i).FIMA-SYS.prevFIMA(i))/(SYS.nowTime(end)-SYS.nowTime(end-1));
              case 'outer'
                diff=(MAT(i).FIMA-SYS.oldFIMA(i))/(SYS.nowTime(end)-SYS.nowTime(end-1));
            end
        end
        if(abs(diff)<OPT.CONV.(loop).value)
          isConverged(i)=true;
        else
          isConverged(i)=false;
        end
      end
    end
end

switch loop
  case 'inner'
    if(all(isConverged))
      SYS.stopInner=true;
      fprintf('%s\n',['** INNER CONVERGENCE ** ' loop ' loop of outer iteration ' ...
        num2str(SYS.ouCntr) ' converged after ' num2str(SYS.inCntr) ' iterations!']);
    end
  case 'outer'
    if(all(isConverged))
      SYS.stopOuter=true;
      fprintf('%s\n',['*** OUTER CONVERGENCE *** ' loop ' loop converged after ' ...
        num2str(SYS.ouCntr) ' iterations!']);
    end
end

return
end

