function neutronBalance(MAT,SYS)
%NEUTRONBALANCE computes neutron balance for the system

matZAI=MAT(1).ZAI;
matNames=MAT(1).nuclideName;
idxAct=isActinide(matZAI);
idxFP=isFP(matZAI);
idxOth=find(~isFP(matZAI)&~isActinide(matZAI));
nMat=length(SYS.IDX.MAT.inFlux);
nNuc=length(matZAI);
nBalance=zeros(nNuc+4,4*nMat);
nCols=5;
for i=1:nMat
    captRate=MAT(SYS.IDX.MAT.inFlux(i)).captRate;
    fissRate=MAT(SYS.IDX.MAT.inFlux(i)).fissRate;
    n2nRate=MAT(SYS.IDX.MAT.inFlux(i)).n2nRate;
    decRate=MAT(SYS.IDX.MAT.inFlux(i)).activity;
    atoms=MAT(SYS.IDX.MAT.inFlux(i)).N(:,end)*1.0E24;
    %Nuclides
    nBalance(1:nNuc,nCols*(i-1)+1)=atoms;
    nBalance(1:nNuc,nCols*(i-1)+2)=captRate;
    nBalance(1:nNuc,nCols*(i-1)+3)=fissRate;
    nBalance(1:nNuc,nCols*(i-1)+4)=n2nRate;
    nBalance(1:nNuc,nCols*(i-1)+5)=decRate;
    %Actinides
    nBalance(nNuc+1,nCols*(i-1)+1)=sum(atoms(idxAct));
    nBalance(nNuc+1,nCols*(i-1)+2)=sum(captRate(idxAct));
    nBalance(nNuc+1,nCols*(i-1)+3)=sum(fissRate(idxAct));
    nBalance(nNuc+1,nCols*(i-1)+4)=sum(n2nRate(idxAct));
    nBalance(nNuc+1,nCols*(i-1)+5)=sum(decRate(idxAct));
    %FPs
    nBalance(nNuc+2,nCols*(i-1)+1)=sum(atoms(idxFP));
    nBalance(nNuc+2,nCols*(i-1)+2)=sum(captRate(idxFP));
    nBalance(nNuc+2,nCols*(i-1)+3)=sum(fissRate(idxFP));
    nBalance(nNuc+2,nCols*(i-1)+4)=sum(n2nRate(idxFP));
    nBalance(nNuc+2,nCols*(i-1)+5)=sum(decRate(idxFP));
    %Others
    nBalance(nNuc+3,nCols*(i-1)+1)=sum(atoms(idxOth));
    nBalance(nNuc+3,nCols*(i-1)+2)=sum(captRate(idxOth));
    nBalance(nNuc+3,nCols*(i-1)+3)=sum(fissRate(idxOth));
    nBalance(nNuc+3,nCols*(i-1)+4)=sum(n2nRate(idxOth));
    nBalance(nNuc+3,nCols*(i-1)+5)=sum(decRate(idxOth));
    %Total
    nBalance(nNuc+4,nCols*(i-1)+1)=sum(atoms);
    nBalance(nNuc+4,nCols*(i-1)+2)=sum(captRate);
    nBalance(nNuc+4,nCols*(i-1)+3)=sum(fissRate);
    nBalance(nNuc+4,nCols*(i-1)+4)=sum(n2nRate);
    nBalance(nNuc+4,nCols*(i-1)+5)=sum(decRate);
end
matZAI(~any(nBalance(:,2:4),2))=[];
matNames(~any(nBalance(:,2:4),2))=[];
nBalance(~any(nBalance(:,2:4),2),:)=[];

idxAct=find(isActinide(matZAI));
idxFP=find(isFP(matZAI));
idxOth=find(~isFP(matZAI)&~isActinide(matZAI));

%[~,B]=sort(sum(nBalance(1:end-4,2:4),2),'DESCEND');
B=vertcat(idxOth,idxAct,idxFP);
[nNuc,nMat]=size(nBalance);
nNuc=nNuc-4;nMat=nMat/nCols;
nBalance=nBalance([B;[nNuc+1:nNuc+4]'],:);
matNames={matNames{B},'Actinides','FPs','Others','Total'};
matZAI=matZAI(B);
matZAI=[matZAI;0;0;0;0];
names={MAT(SYS.IDX.MAT.inFlux).name};

fid=fopen('neutronBalance.txt','w+');
A=repmat({'Atoms','Capture','Fission','(n,2n)','Decay'},1,nMat);
fprintf(fid,['%-10s%-9.5f' repmat(['%-' num2str(13*nCols) 's'],1,nMat) '\n'],'Nu-bar',SYS.RR.NU{2},names{:});
fprintf(fid,'%s\n',repmat('_',1,19+13*nCols*nMat));
fprintf(fid,['%-10s%-9s' repmat(repmat('%-13s',1,nCols),1,nMat) '\n'],'Nuclide','ZAI',A{:});
fprintf(fid,'%s\n',repmat('_',1,19+13*nCols*nMat));
for i=nNuc+1:nNuc+4
    fprintf(fid,['%-10s%-9s' repmat(repmat('%-13G',1,nCols),1,nMat) '\n'],matNames{i},'',nBalance(i,:));
end
fprintf(fid,'%s\n',repmat('_',1,19+13*nCols*nMat));
for i=1:length(idxOth)
    fprintf(fid,['%-10s%-9d' repmat(repmat('%-13G',1,nCols),1,nMat) '\n'],matNames{i},matZAI(i),nBalance(i,:));
end
fprintf(fid,'%s\n',repmat('_',1,19+13*nCols*nMat));
for i=length(idxOth)+1:length(idxOth)+length(idxAct)
    fprintf(fid,['%-10s%-9d' repmat(repmat('%-13G',1,nCols),1,nMat) '\n'],matNames{i},matZAI(i),nBalance(i,:));
end
fprintf(fid,'%s\n',repmat('_',1,19+13*nCols*nMat));
for i=length(idxOth)+length(idxAct)+1:nNuc
    fprintf(fid,['%-10s%-9d' repmat(repmat('%-13G',1,nCols),1,nMat) '\n'],matNames{i},matZAI(i),nBalance(i,:));
end
fclose(fid);

return
end
