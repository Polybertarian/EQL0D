classdef EQL0D
    %EQL0D Case for EQL0D calculation
    
    properties
        name='eql';
        nCores=1;
        verboseMode=false;
        printAndQuit=false;
        restartCalc=false;
        resetCounters=false;
        debugMode=false;
        nuclearDataLibrary='endfb7';
        innerCounter=0;
        outerCounter=0;
        
        OPT
        MAT=[];
        REP
        IDX
        REACT
        REDOX
        
        integralFlux=0.0;
        KEFF=[];
        KINF=[];
    end
    
    methods
        function obj=EQL0D(casename,varargin)%nCores,verbose,print,restart,reset,debug)
            obj.name=casename;
            switch nargin
                case 2
                    obj.nCores=varargin{1};
                case 7
                    obj.nCores=varargin{1};
                    obj.verboseMode=varargin{2};
                    obj.printAndQuit=varargin{3};
                    obj.restartCalc=varargin{4};
                    obj.resetCounters=varargin{5};
                    obj.debugMode=varargin{6};
            end
            run('defaultConfig.m')
            obj.OPT=OPT;
            return
        end
        function obj=set.OPT(obj,OPT)
            obj.OPT=OPT;
            return
        end
        function obj=set.nuclearDataLibrary(obj,library)
            obj.nuclearDataLibrary=library;
        end
        function obj=addMat(obj,name,type,dens,vol,temp,zai,comp)
            load(obj.nuclearDataLibrary);
            obj.MAT(end+1)=Mat(name,type,dens,vol,temp,zai,comp);
            return
        end
        function obj=addRep(obj,name,srcMat,dstMat,elements,share,rate,type)
            obj.REP(end+1)=Rep(name,srcMat,dstMat,elements,share,rate,type);
            return
        end
    end
end

