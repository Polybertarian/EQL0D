function keepMtx = createKeepMatrix(repStream,srcMat,dstMat,totalMtx,burnMat)
%keepMtx = CREATEKEEPMATRIX(MAT,REP) prepares the 'keep'-type reprocessing
%matrices
% if(SYS.debugMode)
%   fprintf('%s\n','*** CONT *** Adding continuous processing streams...');
% end
%
% if(SYS.debugMode)
%   fprintf('%s\n',['*** CONT *** Adding processing stream ' REP.name '...']);
% end
% if(SYS.debugMode)
%   fprintf('%s\n','*** CONT *** Continuous processing streams added!');
% end

if(ismember(repStream.mode,{'keepAFPM','keepAM','keepTotM'}))
  mode='mass';
elseif(ismember(repStream.mode,{'keepAFPA','keepAA','keepTotA'}))
  mode='atoms';
end

dstMatIdx=find(burnMat==repStream.dstMatIdx);
switch repStream.mode
  case {'keepAM','keepAA'}
    nuclToBeRepl=isActinide(dstMat.burnZAI); 
  case {'keepAFPM','keepAFPA'}
    nuclToBeRepl=(isFP(dstMat.burnZAI)|isActinide(dstMat.burnZAI));
  case {'keepTotM','keepTotA'}
    nuclToBeRepl=true(size(dstMat.burnZAI));
end
nuclToBeReplDst=dstMatIdx(nuclToBeRepl);
feedNuc=ismember(dstMat.burnZAI,repStream.elements); % feed nuclides in destination mat
feedNucInDst=dstMatIdx(feedNuc);

% share between feed nuclides
if repStream.srcMatIdx~=0 
  srcMatIdx=find(burnMat==repStream.srcMatIdx);
  feedNuc=ismember(srcMat.burnZAI,repStream.elements);
  feedNucInSrc=srcMatIdx(feedNuc);
  switch mode
    case 'mass'
      share=srcMat.massDens(ismember(srcMat.burnZAI,repStream.elements)); %share taken from source material
    case 'atoms'
      share=srcMat.atDens(ismember(srcMat.burnZAI,repStream.elements));
  end
  if sum(share)==0 
    share=repStream.share;
  end
else
  share=repStream.share;
end
share=share/sum(share);

keepMtx=0.0*totalMtx; % empty matrix of the right size
switch mode
  case 'mass'
    mFeedNuc=repStream.dstNucMass; % masses of feed nuclides
    mReplNuc=dstMat.atomicMass(nuclToBeRepl); % masses of replaced nuclides
    keepMtx(feedNucInDst,nuclToBeRepl)=-share*(sum(totalMtx(nuclToBeReplDst,nuclToBeRepl).*...
      repmat(mReplNuc,1,nnz(nuclToBeRepl)))./sum((mFeedNuc.*share)));
  case 'atoms'
    keepMtx(feedNucInDst,nuclToBeRepl)=-share*sum(totalMtx(nuclToBeReplDst,nuclToBeRepl));
end

if ~isempty(srcMat)
  keepMtx(feedNucInSrc,nuclToBeRepl)=-srcMat.volume/dstMat.volume*keepMtx(feedNucInDst,nuclToBeRepl);
end

return
end

