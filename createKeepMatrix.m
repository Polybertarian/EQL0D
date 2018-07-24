function keepMtx = createKeepMatrix(REP,totalMtx,srcMat,dstMat)
%KeepMtx = CREATEKEEPMATRIX(MAT,REP) prepares the 'keep'-type reprocessing
%matrices
% if(SYS.debugMode)
%   fprintf(SYS.FID.log,'%s\n','*** CONT *** Adding continuous processing streams...');
% end
%
% if(SYS.debugMode)
%   fprintf(SYS.FID.log,'%s\n',['*** CONT *** Adding processing stream ' REP.name '...']);
% end
% if(SYS.debugMode)
%   fprintf(SYS.FID.log,'%s\n','*** CONT *** Continuous processing streams added!');
% end

if(ismember(REP.mode,{'keepAFPM','keepAM','keepTotM'}))
  mode='mass';
elseif(ismember(REP.mode,{'keepAFPA','keepAA','keepTotA'}))
  mode='atoms';
end

switch REP.mode
  case {'keepAM','keepAA'}
    replNuc=isActinide(srcMat.burnZAI);
  case {'keepAFPM','keepAFPA'}
    replNuc=(isFP(srcMat.burnZAI)|isActinide(srcMat.burnZAI));
  case {'keepTotM','keepTotA'}
    replNuc=true(size(srcMat.burnZAI));
end
nucInDst=(matZAI==dst); % nuclides in destination mat
feedNucInDst=ismember(srcMat.burnZAI,REP.elements)&nucInDst; % feed nuclides in destination mat
replNucInDst=replNuc&nucInDst; % replaced nuclides in destination mat

% share between feed nuclides
if(REP.srcMatIdx~=0)
  switch mode
    case 'mass'
      share=srcMat.massDens(ismember(srcMat.ZAI,srcMat.burnZAI(feedNucInDst))); %share taken from source material
    case 'atoms'
      share=srcMat.atDens(ismember(srcMat.ZAI,srcMat.burnZAI(feedNucInDst)));
  end
  if(sum(share)==0)
    share=REP.share;
  end
else
  share=REP.share;
end
share=share/sum(share);

keepMtx=0.0*totalMtx; % empty matrix of the right size
switch mode
  case 'mass'
    mFeedNuc=dstMat.atomicMass(ismember(dstMat.ZAI,srcMat.burnZAI(feedNucInDst))); % masses of feed nuclides
    mReplNuc=dstMat.atomicMass(ismember(dstMat.ZAI,srcMat.burnZAI(replNuc&nucInDst))); % masses of replaced nuclides
    keepMtx(feedNucInDst,replNuc)=-share*(sum(totalMtx(replNucInDst,replNuc).*...
      repmat(mReplNuc,1,nnz(replNuc)))./sum((mFeedNuc.*share)));
  case 'atoms'
    keepMtx(feedNucInDst,replNuc)=-share*sum(totalMtx(replNucInDst,replNuc));
end

return
end

