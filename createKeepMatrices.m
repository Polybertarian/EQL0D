function SYS = createKeepMatrices(MAT,REP,SYS)
%SYS = CREATEKEEPMATRICES(MAT,REP,SYS) prepares the burn-up, decay, and reprocessing matrices for
%the system using the materials vector MAT, reprocessing streams vector REP and system data structure SYS.

%%% Remove old MTX+ create index and vector
SYS.MTX.rep(1,:)=[];
MAT(src).burnZAI=vertcat(SYS.MAT(src).burnZAI{2,:});
matZAI=vertcat(SYS.IDX.matZAI{2,:});
sz=length(MAT(src).burnZAI);

% if(SYS.debugMode)
%   fprintf(SYS.FID.log,'%s\n','*** CONT *** Adding continuous processing streams...');
% end
%
% if(SYS.debugMode)
%   fprintf(SYS.FID.log,'%s\n',['*** CONT *** Adding processing stream ' REP(k).name '...']);
% end

%%% Add continuous processing constants
for k=SYS.IDX.REP.contKeep
  dst=REP(k).dstMatIdx;
  src=REP(k).srcMatIdx;
  tmpMtx=MAT(dst).normBurnMtx; %%% Temporary rep matrix sum
  for i=MAT(dst).streams.cont(MAT(dst).streams.cont<k)
    tmpMtx=tmpMtx+REP(i).repMtx;
  end
  switch REP(k).mode
    case {'keepAM','keepAA'}
      replNuc=isActinide(MAT(src).burnZAI);
    case {'keepAFPM','keepAFPA'}
      replNuc=(isFP(MAT(src).burnZAI)|isActinide(MAT(src).burnZAI));
    case {'keepTotM','keepTotA'}
      replNuc=true(size(MAT(src).burnZAI));
  end
  nucInDst=(matZAI==dst); % nuclides in destination mat
  feedNucInDst=ismember(MAT(src).burnZAI,REP(k).elements)&nucInDst; % feed nuclides in destination mat
  replNucInDst=replNuc&nucInDst; % replaced nuclides in destination mat
  
  % share between feed nuclides
  if(REP(k).srcMatIdx~=0)
    share=MAT(src).massDens(ismember(MAT(src).ZAI,MAT(src).burnZAI(feedNucInDst)));
    if(sum(share)==0)
      share=REP(k).share;
    end
  else
    share=REP(k).share;
  end
  share=share/sum(share);
  
  if(ismember(REP(k).mode,{'keepAFPM','keepAM','keepTotM'}))
    mFeedNuc=MAT(dst).atomicMass(ismember(MAT(dst).ZAI,MAT(src).burnZAI(feedNucInDst))); % masses of feed nuclides
    mReplNuc=MAT(dst).atomicMass(ismember(MAT(dst).ZAI,MAT(src).burnZAI(replNuc&nucInDst))); % masses of replaced nuclides
    repMtx(feedNucInDst,replNuc)=-share*(sum(tmpMtx(replNucInDst,replNuc).*...
      repmat(mReplNuc,1,nnz(replNuc)))./sum((mFeedNuc.*share)));
  elseif(ismember(REP(k).mode,{'keepAFPA','keepAA','keepTotA'}))
    repMtx(feedNucInDst,replNuc)=-share*sum(tmpMtx(replNucInDst,replNuc));
  end
  % if source is a existing material
  if(~isempty(REP(k).srcMatIdx))
    nucInSrc=(matZAI==src);
    feedNucInSrc=ismember(MAT(src).burnZAI,REP(k).elements)&nucInSrc&ismember(MAT(src).burnZAI,MAT(src).burnZAI(feedNucInDst));
    repMtx(feedNucInSrc,replNuc)=-MAT(dst).volume/MAT(src).volume*repMtx(feedNucInDst,replNuc);
  end
  
  SYS.MTX.rep{2,SYS.IDX.REP.cont==k}=sparse(repMtx);
end

if(SYS.debugMode)
  fprintf(SYS.FID.log,'%s\n','*** CONT *** Continuous processing streams added!');
end
return
end

