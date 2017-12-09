classdef Mat < NuclearObject
    %MAT Material for EQL0D
    properties
        name='fuel'
        N;               % Number of atoms per nuclide (1e24 units)
        volume;          % Volume
        intFlux=0.0;
        isBurned=false;  % is the material being depleted?
        isInFlux=false;  % is the material in the neutron flux?
        isStr=false;     % is the material involved in Streams?
        isCont=false;    % is the material continuously processed?
        %%% Cross-sections
        mFissXS;
        mCaptXS;
        mN2nXS;
        mN3nXS;
    end
    properties (SetAccess = private)
        temp;            % Temperature (for Serpent)
        %%% Initial values
        initN;           % Initial number of atoms
        initVolume;      % initial Volume
        initTotMass;        % initial total Mass
        initTotActMass;  % initial total Actinide Mass
        %%% Qualifiers
        streams
    end
    properties (Dependent = true, SetAccess = private)
        %%% Initial values
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
                obj.name=name;
                obj.volume=vol;
                obj.initVolume=vol;
                obj.temp=temp;
                switch type
                    case 1
                        obj.isBurned=true;
                        obj.isInFlux=true;
                        obj.isCont=true;
                    case 2
                        obj.isInFlux=true;
                    case 0
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
                obj.initN=obj.N;
                obj.initTotMass=obj.totMass;
                obj.initTotActMass=obj.totActMass;
                obj.mFissXS=zeros(size(obj.N));
                obj.mCaptXS=zeros(size(obj.N));
                obj.mN2nXS=zeros(size(obj.N));
                obj.mN3nXS=zeros(size(obj.N));
            else
                obj.name='fuel';
                obj.volume=1;
                obj.initVolume=1;
                obj.temp=900;
                obj.N=zeros(size(obj.ZAI));
                obj.initN=obj.N;
            end
        end
        %%% Initial values
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
            totN=sum(obj.N);
        end
        function totActN = get.totActN(obj)
            totActN=sum(obj.N(isActinide(obj.ZAI)));
        end
        function totFPN = get.totFPN(obj)
            totFPN=sum(obj.N(isFP(obj.ZAI)));
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
        %%% Isotope-wise values
        function atDens = get.atDens(obj)
            atDens=obj.N/obj.volume;
        end
        function massDens = get.massDens(obj)
            massDens=1e24*obj.atomicMass.*obj.atDens;
        end
        function actMass = get.actMass(obj)
            actMass=obj.massDens(isActinide(obj.ZAI));
        end
        function molarComp = get.molarComp(obj)
            molarComp=zeros(size(obj.N));
            molarComp(~isElement([8,9,17],obj.ZAI))=obj.N(~isElement([8,9,17],obj.ZAI));
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
            activity=1.0E24*obj.N*log(2)./obj.halfLife;
            activity(isnan(activity)|isinf(activity))=0.0;
        end
        function decayHeat = get.decayHeat(obj)
            decayHeat=obj.activity.*obj.decayEnergy;
        end
        %%% Other
        function denat=get.isDenatured(obj)
            denat=(obj.N(obj.find(922330))+0.6*obj.N(obj.find(922350)))<0.12*sum(obj.N(obj.find(92)));
        end
        function halideExcess = get.halideExcess(obj)
            halideExcess=sum(obj.oxState.*obj.atDens);
        end
        function mainH=get.mainHalide(obj)
            if(obj.atDens(obj.find(90190))>sum(obj.atDens(obj.find([170350 170370]))))
                mainH=obj.find(90190);
            else
                mainH=obj.find([170350 170370]);
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
                case 'fix'
                    fprintf(fid,'%s\n',['mat ' obj.name ' sum burn 1 vol ' num2str(obj.volume,'%.16E') ' fix ' suffix ' ' num2str(obj.temp)]);
                    nucIdx=find(ones(size(obj.ZAI)))';
                case 'nofix'
                    fprintf(fid,'%s\n',['mat ' obj.name ' sum burn 1 vol ' num2str(obj.volume,'%.16E')]);
                    nucIdx=find(hasNuclearData(obj.ZAI))';
                otherwise
                    fprintf(fid,'%s\n',['mat ' obj.name ' sum burn 1 vol ' num2str(obj.volume,'%.16E')]);
                    nucIdx=find(hasNuclearData(obj.ZAI)&obj.N>0.0)';
            end
            for i=nucIdx
                if(hasNuclearData(obj.ZAI(i)))
                    fprintf(fid,'%+11s   %.16E\n',[char(obj.nuclideName(i)) '.' suffix],obj.atDens(i));
                    %elseif(ismember(obj.ZAI(i),[70100,70110,70130,70160,70170,290702,491162,491182,511242,511262,611522,631522,721792]))
                    %    fprintf(fid,'%+11u   %.16E\n',obj.ZAI(i),obj.atDens(i));
                else
                    %fprintf(fid,'%+11s   %.16E\n',char(obj.nuclideName(i)),obj.atDens(i));
                    fprintf(fid,'%+11u   %.16E\n',obj.ZAI(i),obj.atDens(i));
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
            matN=obj.N*1e24;
            matADens=obj.atDens;
            matV=obj.volume;
            matMDens=obj.massDens;
            matM=matMDens*matV/1000.; %kg
            matComp=obj.molarComp;
            matA=obj.activity/1e12;   %TBq
            matH=obj.decayHeat/1e6;        %MW
            matName=obj.nuclideName;
            idxAct=isActinide(matZAI);
            idxFP=isFP(matZAI);
            idxExist=find(isProduced(matZAI));
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
            printFmt.content=['%-9s%-9u' repmat('%-13.6G',1,11) '\n'];
            printFmt.line='%s\n';
            printContent=struct('header',[],'line',[]);
            printContent.header={{'Nuclide','ZAI','Atomic dens.','Tot. atoms',...
                'Molar frac.','Mass dens.','Tot. mass','Activity','Dec. heat',...
                'Toxicity','Capture','Fissions','(n,2n)'},{'','Units','[/b.cm^3]',...
                '[atoms]','[mol%]','[g/cm^3]','[kg]','[TBq]','[MW]','[Sv]','[/s]',...
                '[/s]','[/s]'}};
            printContent.line=repmat('_',1,15*11);
            fprintf(fid,printFmt.line,printContent.line);
            fprintf(fid,printFmt.header,printContent.header{1}{:});
            fprintf(fid,printFmt.header,printContent.header{2}{:});
            fprintf(fid,printFmt.line,printContent.line);
            fprintf(fid,['%-18s' repmat('%-13.6G',1,11) '\n'],'FPs',sum(matADens(idxFP)),...
                sum(matN(idxFP)) ,sum(matComp(idxFP)) ,sum(matMDens(idxFP)) ,sum(matM(idxFP)),...
                sum(matA(idxFP)) ,sum(matH(idxFP)) ,sum(matTox(idxFP)) ,sum(matCapt(idxFP)),...
                sum(matFiss(idxFP)) , sum(matN2n(idxFP)) );
            fprintf(fid,['%-18s' repmat('%-13.6G',1,11) '\n'],'Actinides',sum(matADens(idxAct)),...
                sum(matN(idxAct)),sum(matComp(idxAct)),sum(matMDens(idxAct)),sum(matM(idxAct)),...
                sum(matA(idxAct)),sum(matH(idxAct)),sum(matTox(idxAct)),sum(matCapt(idxAct)),...
                sum(matFiss(idxAct)),sum(matN2n(idxAct)));
            fprintf(fid,['%-18s' repmat('%-13.6G',1,11) '\n'],'Total',sum(matADens),...
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
        function idx=find(obj,ZAI)
            if(all(ZAI<1000))
                idx=find(isElement(ZAI,obj.ZAI));
            else
                idx=find(ismember(obj.ZAI,ZAI));
            end
        end
        function mFrac=mFrac(obj,IDX)
            mFrac=obj.massDens(IDX);
            mFrac=mFrac/sum(mFrac);
        end
        function aFrac=aFrac(obj,IDX)
            aFrac=obj.atDens(IDX);
            aFrac=aFrac/sum(aFrac);
        end
        function avMass=avMass(obj,IDX)
           avMass=sum(obj.aFrac(IDX).*obj.atomicMass(IDX));
        end
    end
end


