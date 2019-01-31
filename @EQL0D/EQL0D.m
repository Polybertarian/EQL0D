classdef EQL0D
%EQL0D Case for EQL0D calculation

properties
caseName;
iterMode;
nCores=1;
verboseMode=false;
printAndQuit=false;
restartCalc=false;
resetCounters=false;
nuclearDataLibrary;
debugMode=false;
innerCounter=0;
outerCounter=0;

OPT
MAT=[];
REP=[];
IDX
REACT
REDOX
DAT
REA

integralFlux=0.0;
KEFF=[];
KINF=[];
end

methods
function obj=EQL0D(library,iterMode,nCycles,nSteps,cycleLength,printOptions)%nCores,verbose,print,restart,reset,debug)
    obj.nuclearDataLibrary=library
    obj.DAT=load(obj.nuclearDataLibrary,'DAT')
end
function obj=set.OPT(obj,OPT)
    obj.OPT=OPT;
end
function obj=addMat(obj,name,type,dens,vol,temp,zai,comp)
    global DAT=obj.DAT.DAT;
    obj.MAT(end+1)=Mat(name,type,dens,vol,temp,zai,comp);
end
function obj=addRep(obj,name,srcMat,dstMat,elements,share,rate,type)
    obj.REP(end+1)=Rep(name,srcMat,dstMat,elements,share,rate,type);
end
function obj=setOptions(obj,caseName,nCores,verboseMode,printAndQuit,restartCalc,resetCounters,debugMode)
    obj.caseName=caseName;
    obj.nCores=nCores;
    obj.verboseMode=verboseMode;
    obj.printAndQuit=printAndQuit;
    obj.restartCalc=restartCalc;
    obj.resetCounters=resetCounters;
    obj.debugMode=debugMode;
end
end
end
