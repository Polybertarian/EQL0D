%%% Default values
OPT.serpentPath='nice -n 10 runsss /afs/psi.ch/project/fast_lrs/workspace/COD/SERP/Serpent2/EQL0D/2.1.26/sss2';
OPT.cycleLength=365;
OPT.nCycles=9;
OPT.nSteps=[50 2000];
OPT.iterMode='equilibrium';
OPT.printSteps=false;
OPT.printStepsBB=false;
OPT.printCycles=true;
OPT.defaultDataLibrary='endfb7';
OPT.PCC=false;

%%% Reactivity control
OPT.reactControl=false;
OPT.REA.mode='replace';
OPT.REA.allowNegative=true;
OPT.REA.targetKeff=1.0;
OPT.REA.upNuclides=922330;
OPT.REA.upFraction=[];
OPT.REA.downNuclides=922330;
OPT.REA.downFraction=[];
OPT.REA.replNuclides=902320;
OPT.REA.replFraction=[];
OPT.REA.feedMat='feed';
OPT.REA.targetMat='fuel';
OPT.REA.tol=300; %pcm
OPT.REA.maxIter=20; 

%%% Redox control
OPT.redoxControl=false;
OPT.REDOX.materials='fuel';
OPT.REDOX.replaceMode='replaceMass';
OPT.REDOX.replaceWith=922380;

%%% Misc
OPT.writeMail=true;
OPT.matWriteStyle='nofix';
OPT.keepFiles=true;

%%% Convergence
OPT.CONV.inner.criterion='maxNucDatRelDiff';
OPT.CONV.inner.value=1e-3;
OPT.CONV.inner.cutoff=1e-10;

OPT.CONV.outer.criterion='maxNucDatRelDiff';
OPT.CONV.outer.value=1e-3;
OPT.CONV.outer.cutoff=1e-10;