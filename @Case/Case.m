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
    
    integralFlux=0.0;
    KEFF=[];
    KINF=[];
  end
  
  methods
    function obj=Case(name,verbose,print)
      obj.name=name;
      obj.verboseMode=verbose;
      obj.printAndQuit=print;
      
      return
    end
  end
end

