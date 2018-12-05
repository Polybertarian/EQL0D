classdef Mat
    %MAT Material for EQL0D
    properties
        name
        N                % Number of atoms per nuclide (1e24 units)
        volume
        intFlux=0.0
        mFissXS          % Cross-sections
        mCaptXS
        mN2nXS
        mN3nXS
        burnMtx={[],[],[]};
        burnIdx
        streams=struct('cont',[],'batch',[]);
    end
    properties (SetAccess = immutable)
        atomicMass
        decayEnergy
        halfLife
        ZAI
        oxState
        ingToxicity
        inhToxicity
        hasNucData
        nuclideName
        defDecMtx
        temp;            % Temperature (for Serpent)
        isBurned=false;  % is the material being depleted?
        isInFlux=false;  % is the material in the neutron flux?
        isStr=false;     % is the material involved in Streams?
        isCont=false;    % is the material continuously processed?
    end
    properties (Dependent = true, SetAccess = private)
        %%% Initial values
        initN;           % Initial number of atoms
        initVolume;      % initial Volume
        initTotMass;     % initial total Mass
        initTotActMass;  % initial total Actinide Mass
        initActN         % initial actinide atoms
        initTotActN      % initial total actinide atoms
        initTotN
        %%% Values per isotope
        atDens           % atomic densities
        massDens         % mass densities
        molarComp        % molar composition
        actMass          % actinide masses
        FPMass           % FP masses
        ingTox           % ingestion toxicity
        activity         % activity
        decayHeat        % decay heat
        %%% Reaction rates
        fissRate         % fission rate
        captRate         % capture rate
        n2nRate          % (n,2n) rate
        n3nRate          % (n,3n) rate
        %%% Total values
        totN             % total number of atoms
        totFPN           % total number of FP atoms
        totActN          % total number of Actinide atoms
        totFPMass        % total FPs mass
        totActMass       % total Actinide mass
        totMass          % total Mass
        %%%%
        burnZAI
        decMtx
        %%%% Other
        FIMA             % burnup (% FIMA)
        FPFrac           % burnup (FP/FP+Ac) mass ratio
        isDenatured
        mainHalide
        halideExcess
    end
    methods
        %Constructor
        function obj = Mat(name,type,dens,vol,temp,zai,comp)
            global DAT
            if exist('DAT','var')
                obj.burnIdx=isProduced(DAT.libraryName,DAT.ZAI0);
                obj.ZAI=DAT.ZAI0(obj.burnIdx);
                obj.atomicMass=DAT.AMASS(obj.burnIdx);
                obj.decayEnergy=DAT.Q(obj.burnIdx);
                obj.halfLife=DAT.T12(obj.burnIdx);
                obj.oxState=valenceStates(obj.ZAI);
                [obj.ingToxicity,obj.inhToxicity]=ingAndInhTox(obj.ZAI);
                obj.hasNucData=hasXSData(obj.ZAI);
                obj.nuclideName=ZAI2Name(obj.ZAI);
                defDecMtx=sparse(DAT.decayMatrix);
                defDecMtx=defDecMtx(:,obj.burnIdx);
                defDecMtx=defDecMtx(obj.burnIdx,:);
                obj.defDecMtx=defDecMtx;
                obj.burnIdx=ismember(obj.ZAI,obj.ZAI);
            else
                error('Nuclear data library not found.')
            end
            if(nargin > 0)
                if(~iscolumn(zai))
                    zai=zai';
                end
                if(~iscolumn(comp))
                    comp=comp';
                end
                if(length(zai)~=length(comp))
                    error('Nuclide vector and composition vector length mismatch!')
                end
                idx=find(~ismember(zai,obj.ZAI));
                if(~isempty(idx))
                    error(['Nuclide(s) ' num2str(zai(idx)) ' not found in library.'])
                end
                [zai,idx]=sort(zai);
                comp=comp(idx);
                obj.name=name;
                obj.volume=vol;
                obj.temp=temp;
                if(~isstr(type))
                    if(type==0)
                        type='decay';
                    elseif(type==1)
                        type='burned';
                    elseif(type==2)
                        type='inFlux';
                    end
                end
                switch type
                    case 'burned'
                        disp('lol')
                        obj.isBurned=true;
                        obj.isInFlux=true;
                        obj.isCont=true;
                    case 'inFlux'
                        obj.isInFlux=true;
                    case 'decay'
                        obj.isCont=true;
                        obj.isStr=true;
                end
                obj.N=zeros(size(obj.ZAI));
                if(dens==0)%'sum' in serpent
                    if(~all(comp<0))
                        obj.N(ismember(obj.ZAI,zai))=obj.volume*comp;
                    elseif(all(comp)==0||isempty(comp))
                        obj.N=zeros(size(obj.ZAI));
                    else
                        error('Error: Negative composition given with 0 density')
                    end
                elseif(dens<0) %mass density given
                    if(all(comp>0))
                        comp=comp/sum(comp);
                        obj.N(ismember(obj.ZAI,zai))=obj.volume*-dens*comp/...
                            sum(1e24*obj.atomicMass(ismember(obj.ZAI,zai)).*comp);
                    elseif(all(comp<0))
                        comp=comp/sum(comp);
                        obj.N(ismember(obj.ZAI,zai))=comp*-dens*obj.volume/...
                            sum(1e24*obj.atomicMass(ismember(obj.ZAI,zai)));
                    else
                        error('Error: Mixed positive and negative values in composition!')
                    end
                else %atomic density given
                    if(all(comp>0))
                        comp=comp/sum(comp);
                        obj.N(ismember(obj.ZAI,zai))=comp*dens*obj.volume;
                    elseif(all(comp<0))
                        comp=comp/sum(comp);
                        tmp = comp./(1e24*obj.atomicMass);
                        obj.N(ismember(obj.ZAI,zai))=obj.volume*comp*dens.*tmp/sum(tmp);
                    else
                        error('Error: Mixed positive and negative values in composition!')
                    end
                end
                obj.mFissXS=zeros(size(obj.N));
                obj.mCaptXS=zeros(size(obj.N));
                obj.mN2nXS=zeros(size(obj.N));
                obj.mN3nXS=zeros(size(obj.N));
            else
                obj.name='fuel';
                obj.volume=1;
                obj.temp=900;
                obj.oldN=[];
                obj.N=zeros(size(obj.ZAI));
            end
        end
        %%% Initial values
        function initN = get.initN(obj)
            initN=obj.N(:,1);
        end
        function initVolume = get.initVolume(obj)
            initVolume=obj.volume(1);
        end
        function initTotMass=get.initTotMass(obj)
            initTotMass=sum(obj.atomicMass.*obj.initN*1E24);
        end
        function initTotActMass=get.initTotActMass(obj)
            initTotActMass=sum(obj.atomicMass(isActinide(obj.ZAI)).*obj.initActN*1E24);
        end
        function initTotN = get.initTotN(obj)
            initTotN=sum(obj.initN);
        end
        function initActN = get.initActN(obj)
            initActN=obj.initN(isActinide(obj.ZAI));
        end
        function initTotActN = get.initTotActN(obj)
            initTotActN=sum(obj.initActN);
        end
        %%% Total values
        function totN = get.totN(obj)
            totN=sum(obj.N(:,end));
        end
        function totActN = get.totActN(obj)
            totActN=sum(obj.N(isActinide(obj.ZAI),end));
        end
        function totFPN = get.totFPN(obj)
            totFPN=sum(obj.N(isFP(obj.ZAI),end));
        end
        function totActMass = get.totActMass(obj)
            totActMass=obj.volume*sum(obj.actMass);
        end
        function totMass = get.totMass(obj)
            totMass=obj.volume*sum(obj.massDens);
        end
        function totFPMass=get.totFPMass(obj)
            totFPMass=obj.volume*sum(obj.FPMass);
        end
        function atDens = get.atDens(obj)
            atDens=obj.N(:,end)/obj.volume;
        end
        function massDens = get.massDens(obj)
            massDens=1E24*obj.atomicMass.*obj.atDens;
        end
        function actMass = get.actMass(obj)
            actMass=obj.massDens(isActinide(obj.ZAI));
        end
        function molarComp = get.molarComp(obj)
            molarComp=zeros(size(obj.N(:,end)));
            molarComp(~isElement([8,9,17],obj.ZAI))=obj.N(~isElement([8,9,17],obj.ZAI),end);
            molarComp=100*molarComp/sum(molarComp);
            molarComp(isnan(molarComp))=0.0;
        end
        %%% Reaction rates
        function fissRate=get.fissRate(obj)
            fissRate=obj.mFissXS.*obj.atDens.*obj.intFlux;
        end
        function captRate=get.captRate(obj)
            captRate=obj.mCaptXS.*obj.atDens.*obj.intFlux;
        end
        function n2nRate =get.n2nRate(obj)
            n2nRate=obj.mN2nXS.*obj.atDens.*obj.intFlux;
        end
        function n3nRate =get.n3nRate(obj)
            n3nRate=obj.mN3nXS.*obj.atDens.*obj.intFlux;
        end
        %%% Isotope values
        function ingTox=get.ingTox(obj)
            ingTox=obj.activity.*obj.ingToxicity;
        end
        function activity = get.activity(obj)
            activity=1.0E24*obj.N(:,end)*log(2)./obj.halfLife;
            activity(isnan(activity)|isinf(activity))=0.0;
        end
        function decayHeat = get.decayHeat(obj)
            decayHeat=obj.activity.*obj.decayEnergy;
        end
        %%% Other
        function denat=get.isDenatured(obj)
            denat=(obj.N(obj.find(922330),end)+0.6*obj.N(obj.find(922350),end))<0.12*sum(obj.N(obj.find(92),end));
        end
        function halideExcess = get.halideExcess(obj)
            halideExcess=sum(obj.oxState.*obj.N(:,end));
        end
        function mainH=get.mainHalide(obj)
            idxF=obj.find(90190);
            idxCl=obj.find([170350 170370]);
            if(obj.atDens(idxF)>sum(obj.atDens(idxCl)))
                mainH=idxF;
            else
                if(obj.atDens(idxCl(1))==0)
                    mainH=idxCl(2);
                elseif(obj.atDens(idxCl(2))==0)
                    mainH=idxCl(1);
                else
                    mainH=idxCl;
                end
            end
        end
        function FPMass=get.FPMass(obj)
            FPMass=obj.massDens(isFP(obj.ZAI));
        end
        function FIMA=get.FIMA(obj)
            FIMA=100.0*(obj.initTotActMass-obj.totActMass)/obj.initTotActMass;
        end
        function FPFrac=get.FPFrac(obj)
            FPFrac=100.0*obj.totFPMass/(obj.totFPMass+obj.totActMass);
        end
        %%% Write/Print
        function status = write(obj,matWriteStyle)
            fileName=[obj.name '.serp'];
            %Check if existing file & backup
            if(exist(['./' fileName],'file')==2)
                movefile(fileName,[fileName '.bak'],'f')
            end
            %Open file
            [fid,errmsg]=fopen(['./' fileName],'w');
            if(~isempty(errmsg))
                error(errmsg);
            end
            suffix=[num2str(obj.temp/100,'%02d') 'c'];
            switch matWriteStyle
                case 'dec'
                    fprintf(fid,'%s\n',['mat ' obj.name ' sum burn 1 vol ' num2str(obj.volume,'%.16E') ' fix ' ...
                        obj.name ' 300']);
                    nucIdx=find(~ismember(obj.ZAI,[1:1:111]*1E4)&obj.N(:,end)>0.0)'; %avoid X-nat nuclides
                    suffix='';
                case 'fix'
                    fprintf(fid,'%s\n',['mat ' obj.name ' sum burn 1 vol ' num2str(obj.volume,'%.16E') ' fix ' ...
                        suffix ' ' num2str(obj.temp)]);
                    nucIdx=find(ones(size(obj.ZAI)))';
                    suffix=['.' num2str(obj.temp/100,'%02d') 'c'];
                case 'nofix'
                    fprintf(fid,'%s\n',['mat ' obj.name ' sum burn 1 vol ' num2str(obj.volume,'%.16E')]);
                    nucIdx=find(hasXSData(obj.ZAI))';
                    suffix=['.' num2str(obj.temp/100,'%02d') 'c'];
                otherwise
                    fprintf(fid,'%s\n',['mat ' obj.name ' sum burn 1 vol ' num2str(obj.volume,'%.16E')]);
                    nucIdx=find(hasXSData(obj.ZAI)&obj.N(:,end)>0.0)';
                    suffix=['.' num2str(obj.temp/100,'%02d') 'c'];
            end
            for i=nucIdx
                if(strcmp(matWriteStyle,'dec')||~hasSymbol(obj.ZAI(i)))
                    fprintf(fid,'%-13u%.16E\n',obj.ZAI(i),obj.atDens(i));
                else
                    fprintf(fid,'%-13s%.16E\n',[char(obj.nuclideName(i)) suffix],obj.atDens(i));
                end
            end
            status=fclose(fid);
        end
        function status = printMaterial(obj,SYS,date)
            switch date
                case 'BB'
                    suffix=['c' num2str(SYS.ouCntr,'%03d') '_s' num2str(SYS.inCntr+1,'%03d') '_BB'];
                case 'AB'
                    suffix=['c' num2str(SYS.ouCntr,'%03d') '_s' num2str(SYS.inCntr,'%03d') '_AB'];
                case 'EoC'
                    suffix=['c' num2str(SYS.ouCntr,'%03d') '_AB'];
                otherwise
                    suffix=date;
            end
            [fid,errmsg]=fopen([obj.name '_' suffix '.txt'],'wt'); %%% Open file
            if(~isempty(errmsg))
                error(errmsg)
            end
            matZAI=obj.ZAI;
            matN=obj.N(:,end)*1E24;
            matADens=obj.atDens;
            matV=obj.volume;
            matMDens=obj.massDens;
            matM=matMDens*matV/1000; %kg
            matComp=obj.molarComp;
            matA=obj.activity;   %Bq
            matH=obj.decayHeat/1E6;        %MW
            matName=obj.nuclideName;
            idxAct=isActinide(matZAI);
            idxFP=isFP(matZAI);
            idxExist=find(obj.burnIdx);
            idxAFP=idxFP|idxAct;
            matFiss=obj.fissRate;
            matCapt=obj.captRate;
            matN2n=obj.n2nRate;
            matTox=obj.ingTox;
            if(obj.isDenatured)
                denat='U-Vector denatured.';
            else
                denat='U-Vector NOT denatured.';
            end
            %%% Header
            if(obj.isBurned)
                fprintf(fid,'%s\n',['Material ''' obj.name ''', volume ' num2str(obj.volume,'%.4G')...
                    ' cm3, time ' num2str(SYS.nowTime(end),'%.4G') ' EFPD, burn-up '...
                    num2str(obj.FIMA,'%.3f') ' %FIMA or ' num2str(obj.FPFrac,'%.3f') ' %FP/(FP+A). ' denat]);
            else
                fprintf(fid,'%s\n',['Material ''' obj.name ''', volume ' num2str(obj.volume,'%.4G')...
                    ' cm3, time ' num2str(SYS.nowTime(end),'%.4G') ' EFPD.']);
            end
            
            printFmt=struct('header',[],'content',[],'line',[]);
            printFmt.header=['%-9s%-9s' repmat('%-13s',1,11) '\n'];
            printFmt.content=['%-9s%-9u' repmat('%-13.5E',1,11) '\n'];
            printFmt.line='%s\n';
            printContent=struct('header',[],'line',[]);
            printContent.header={{'Nuclide','ZAI','Atomic dens.','Tot. atoms','Molar frac.','Mass dens.',...
                'Tot. mass','Activity','Dec. heat','Toxicity','Capture','Fissions','(n,2n)'},{'','Units',...
                '[/b.cm^3]','[atoms]','[mol%]','[g/cm^3]','[kg]','[Bq]','[MW]','[Sv]','[/s]','[/s]','[/s]'}};
            printContent.line=repmat('_',1,15*11);
            fprintf(fid,printFmt.line,printContent.line);
            fprintf(fid,printFmt.header,printContent.header{1}{:});
            fprintf(fid,printFmt.header,printContent.header{2}{:});
            fprintf(fid,printFmt.line,printContent.line);
            fprintf(fid,['%-10s%-8s' repmat('%-13.5E',1,11) '\n'],'FPs','20<Z<73',sum(matADens(idxFP)),...
                sum(matN(idxFP)) ,sum(matComp(idxFP)) ,sum(matMDens(idxFP)) ,sum(matM(idxFP)),...
                sum(matA(idxFP)) ,sum(matH(idxFP)) ,sum(matTox(idxFP)) ,sum(matCapt(idxFP)),...
                sum(matFiss(idxFP)) , sum(matN2n(idxFP)) );
            fprintf(fid,['%-10s%-8s' repmat('%-13.5E',1,11) '\n'],'Actinides','Z>88',sum(matADens(idxAct)),...
                sum(matN(idxAct)),sum(matComp(idxAct)),sum(matMDens(idxAct)),sum(matM(idxAct)),...
                sum(matA(idxAct)),sum(matH(idxAct)),sum(matTox(idxAct)),sum(matCapt(idxAct)),...
                sum(matFiss(idxAct)),sum(matN2n(idxAct)));
            fprintf(fid,['%-10s%-8s' repmat('%-13.5E',1,11) '\n'],'A+FPs','Z>20',sum(matADens(idxAFP)),...
                sum(matN(idxAFP)),sum(matComp(idxAFP)),sum(matMDens(idxAFP)),sum(matM(idxAFP)),...
                sum(matA(idxAFP)),sum(matH(idxAFP)),sum(matTox(idxAFP)),sum(matCapt(idxAFP)),...
                sum(matFiss(idxAFP)),sum(matN2n(idxAFP)));
            fprintf(fid,['%-10s%-8s' repmat('%-13.5E',1,11) '\n'],'Total','All',sum(matADens),...
                sum(matN),sum(matComp),sum(matMDens),sum(matM),...
                sum(matA),sum(matH),sum(matTox),sum(matCapt),...
                sum(matFiss),sum(matN2n));
            fprintf(fid,printFmt.line,printContent.line);
            
            %%% Nuclide data
            for k=idxExist'
                fprintf(fid,printFmt.content,char(matName{k}),matZAI(k),matADens(k),matN(k),...
                    matComp(k),matMDens(k),matM(k),matA(k),matH(k),matTox(k),matCapt(k),...
                    matFiss(k),matN2n(k));
            end
            fprintf(fid,printFmt.line,printContent.line);
            for j=1:111
                elemIdx=isElement(j,obj.ZAI);
                elemN=sum(matN(elemIdx));
                elemADens=sum(matADens(elemIdx));
                elemMDens=sum(matMDens(elemIdx));
                elemM=sum(matM(elemIdx));
                elemComp=sum(matComp(elemIdx));
                elemA=sum(matA(elemIdx));
                elemH=sum(matH(elemIdx));
                elemName=ZAI2Name(j);
                elemFiss=sum(matFiss(elemIdx));
                elemCapt=sum(matCapt(elemIdx));
                elemN2n=sum(matN2n(elemIdx));
                elemTox=sum(matTox(elemIdx));
                fprintf(fid,printFmt.content,char(elemName{1}),j,elemADens,elemN,elemComp,...
                    elemMDens,elemM,elemA,elemH,elemTox,elemCapt,elemFiss,elemN2n);
            end
            %%% Close file
            status=fclose(fid);
        end
        %%% Methods
        function idx = find(obj,ZAI)
            if(all(ZAI<1000))
                idx=find(isElement(ZAI,obj.ZAI));
            else
                idx=find(ismember(obj.ZAI,ZAI));
            end
        end
        function mFrac = mFrac(obj,IDX)
            mFrac=obj.massDens(IDX);
            mFrac=mFrac/sum(mFrac);
        end
        function aFrac = aFrac(obj,IDX)
            aFrac=obj.atDens(IDX);
            aFrac=aFrac/sum(aFrac);
        end
        function avMass = avMass(obj,IDX)
            avMass=sum(obj.aFrac(IDX).*obj.atomicMass(IDX));
        end
        function [obj,sol] = redoxControl(obj,halideIdx,nucIdx,mode)
            excess=sum(obj.oxState.*obj.N(:,end));
            if(excess~=0)
                alpha_halide=obj.N(halideIdx,end)/sum(obj.N(halideIdx,end));
                alpha_target=obj.N(nucIdx,end)/sum(obj.N(nucIdx,end));
                switch mode
                    case 'replaceMass'
                        sol=[alpha_halide'*obj.oxState(halideIdx),alpha_target'*obj.oxState(nucIdx);...
                            -obj.avMass(halideIdx),obj.avMass(nucIdx)]\[excess;0];
                    case 'replaceAtom'
                        sol=[alpha_halide'*obj.oxState(halideIdx),alpha_target'*obj.oxState(nucIdx);1,1]\[excess;0];
                    case 'remove'
                        sol=[excess;0];
                end
                obj.N(halideIdx,end)=obj.N(halideIdx,end)+sol(1).*alpha_halide;
                obj.N(nucIdx,end)   =obj.N(nucIdx,end)+sol(2).*alpha_target;
            end
        end
        function dN=lastNChange(obj)
            dN=obj.N(:,end)-obj.N(:,end-1);
        end
        function nBMtx = normBurnMtx(obj,coeffs)
            if(obj.isInFlux)
                nBMtx=obj.burnMtx{3}*coeffs(3);
                for i=1:2
                    if ~isempty(obj.burnMtx{i})
                        nBMtx=nBMtx+obj.burnMtx{i}*coeffs(i);
                    end
                end
                nBMtx = nBMtx*obj.intFlux;
            else
                nBMtx=0.0*obj.decMtx;
            end
        end
        function burnZAI=get.burnZAI(obj)
            burnZAI=obj.ZAI(obj.burnIdx);
        end
        function obj=loadBurnMatrix(obj)
            run(strtrim(ls(['depmtx_' obj.name '*0.m'])));
            clearvars N0 N1 t
            A=sparse(A);
            idx=find(ismember(ZAI,[-1 10 666]));
            A(:,idx)=[]; A(idx,:)=[]; ZAI(idx)=[];
            obj.burnIdx=ismember(obj.ZAI,ZAI);
            obj.burnMtx(1)=[];
            obj.burnMtx{3}=(A-obj.defDecMtx)/obj.intFlux;
        end
        function decMtx=get.decMtx(obj)
            decMtx=obj.defDecMtx;
            decMtx(:,~obj.burnIdx)=[];
            decMtx(~obj.burnIdx,:)=[];
        end
    end
  end
end