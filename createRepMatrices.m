function SYS = createRepMatrices(MAT,REP,SYS)
%SYS = CREATEREPMATRICES(MAT,REP,SYS) prepares the burn-up, decay, and reprocessing matrices for
%the system using the materials vector MAT, reprocessing streams vector REP and system data structure SYS.

%%% Remove old MTX+ create index and vector
SYS.MTX.rep(1,:)=[];
burnZAI=vertcat(SYS.burnZAI{2,:});
matZAI=vertcat(SYS.IDX.matZAI{2,:});
skipped=[]; % skip "keepX" streams
sz=length(burnZAI);

if(SYS.debugMode)
    fprintf(SYS.FID.log,'%s\n','*** CONT *** Adding continuous processing streams...');
end

%%% Add continuous processing constants
for k=SYS.IDX.contStr
    if(SYS.debugMode)
        fprintf(SYS.FID.log,'%s\n',['*** CONT *** Adding processing stream ' REP(k).name '...']);
    end
    if(~REP(k).isKeep)
        if(isempty(REP(k).srcMatIdx)) %%% void source = find destination elements
            nucInDst=(matZAI==REP(k).dstMatIdx);
            replNuc=[]; % vector isotopes in source material
            feedNucInDst=find(ismember(burnZAI,REP(k).elements)&nucInDst);
        elseif(isempty(REP(k).dstMatIdx)) %%% void destination = find source elements
            nucInSrcMat=(matZAI==REP(k).srcMatIdx);
            replNuc=find(ismember(burnZAI,REP(k).elements)&nucInSrcMat);
            feedNucInDst=[];
        else %%% non-void source & destination = find both
            nucInDst=(matZAI==REP(k).dstMatIdx);
            nucInSrcMat=(matZAI==REP(k).srcMatIdx);
            replNuc=find(ismember(burnZAI,REP(k).elements)&nucInSrcMat);
            feedNucInDst=find(ismember(burnZAI,burnZAI(replNuc))&nucInDst);
            if(length(replNuc)>length(feedNucInDst))
                replNuc=find(ismember(burnZAI,burnZAI(feedNucInDst))&nucInSrcMat);
            end
        end
        SYS.MTX.rep{2,SYS.IDX.contStr==k}=spalloc(sz,sz,length(REP(k).elements));
        %%% modify matrix for found elements
        for l=1:length(replNuc)
            if(REP(k).srcMatIdx~=0)
                SYS.MTX.rep{2,SYS.IDX.contStr==k}(replNuc(l),replNuc(l))=-REP(k).rate;
                if(REP(k).dstMatIdx~=0)
                    SYS.MTX.rep{2,SYS.IDX.contStr==k}(feedNucInDst(l),replNuc(l))=...
                        +REP(k).rate*MAT(REP(k).srcMatIdx).volume/MAT(REP(k).dstMatIdx).volume;
                end
            else
                SYS.MTX.rep{2,SYS.IDX.contStr==k}(feedNucInDst(l),feedNucInDst(l))=+REP(k).rate;
            end
        end
    else
        skipped(end+1)=k;
    end
    clearvars l
end

%%% KeepX modes at the end
for k=skipped
    dst=REP(k).dstMatIdx;
    src=REP(k).srcMatIdx;
    tmpMtx=blkdiag(SYS.MTX.burn{2,:}); %%% Temporary rep matrix sum
    for i=SYS.IDX.contStr(SYS.IDX.contStr<k)
        tmpMtx=tmpMtx+SYS.MTX.rep{2,i};
    end
    switch REP(k).mode
        case {'keepAM','keepAA'}
            replNuc=isActinide(burnZAI);
        case {'keepAFPM','keepAFPA'}
            replNuc=(isFP(burnZAI)|isActinide(burnZAI));
        case {'keepTotM','keepTotA'}
            replNuc=true(size(burnZAI));
    end
    nucInDst=(matZAI==dst); % nuclides in destination mat
    feedNucInDst=ismember(burnZAI,REP(k).elements)&nucInDst; % feed nuclides in destination mat
    replNucInDst=replNuc&nucInDst; % replaced nuclides in destination mat
    
    % allocate empty matrix
    repMtx=zeros(sz,sz);
    
    % share between feed nuclides
    if(REP(k).srcMatIdx~=0)
        share=MAT(src).massDens(ismember(MAT(src).ZAI,burnZAI(feedNucInDst)));
        if(sum(share)==0)
            share=REP(k).share;
        end
    else
        share=REP(k).share;
    end
    share=share/sum(share);
    
    if(ismember(REP(k).mode,{'keepAFPM','keepAM','keepTotM'}))
        mFeedNuc=MAT(dst).atomicMass(ismember(MAT(dst).ZAI,burnZAI(feedNucInDst))); % masses of feed nuclides
        mReplNuc=MAT(dst).atomicMass(ismember(MAT(dst).ZAI,burnZAI(replNuc&nucInDst))); % masses of replaced nuclides
        repMtx(feedNucInDst,replNuc)=-share*(sum(tmpMtx(replNucInDst,replNuc).*...
            repmat(mReplNuc,1,nnz(replNuc)))./sum((mFeedNuc.*share)));
    elseif(ismember(REP(k).mode,{'keepAFPA','keepAA','keepTotA'}))
        repMtx(feedNucInDst,replNuc)=-share*sum(tmpMtx(replNucInDst,replNuc));
    end
    % if source is a existing material
    if(~isempty(REP(k).srcMatIdx))
        nucInSrc=(matZAI==src);
        feedNucInSrc=ismember(burnZAI,REP(k).elements)&nucInSrc&ismember(burnZAI,burnZAI(feedNucInDst));
        repMtx(feedNucInSrc,replNuc)=-MAT(dst).volume/MAT(src).volume*repMtx(feedNucInDst,replNuc);
    end
    
    SYS.MTX.rep{2,SYS.IDX.contStr==k}=sparse(repMtx);
end

if(SYS.debugMode)
    fprintf(SYS.FID.log,'%s\n','*** CONT *** Continuous processing streams added!');
end
return
end

