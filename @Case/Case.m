classdef Case
  %CASE Case for EQL0D calculation
  
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
    MAT
    IDX
    REACT
    REDOX
    
    integralFlux=0.0;
    KEFF=[];
    KINF=[];
  end
  
  methods
    function obj=Case(name,nCores,verbose,print,restart,reset,debug)
      obj.name=name;
      obj.nCores=nCores;
      obj.verboseMode=verbose;
      obj.printAndQuit=print;
      obj.restartCalc=restart;
      obj.resetCounters=reset;
      obj.debugMode=debug;
      return
    end
    function obj=set.OPT(obj,OPT)
      obj.OPT=OPT;
    end
  end
end

