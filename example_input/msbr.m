load('endfb7')
OPT.nSteps=6; OPT.nCycles=30;
OPT.cycleLength=[15 20 25 32 40 51 65 82 104 132 168 213 270 343 435 552 700 888 1127 1430 1815 2303 2923 3709 4706 5972 7578 9616 12203 15485];
OPT.iterMode='steps';
OPT.reactControl=true;
OPT.PCC=true;

%Reactivity control
OPT.REA.feedMat='feed';
OPT.REA.mode='addMass';
OPT.REA.upNuclides=92;
OPT.REA.downNuclides=[922320 922330 922340 922350];
OPT.REA.downFraction=[0.0002406 0.9921 0.0076598 1.6178E-06];

ThVec=902320;ThF=1;
initfuel=[68.5 31.3 0.2]; UShare=0.176; UVec=[922330]; UFrac=[1];
fuelcomp=[initfuel(1:2) UShare]; fuelcomp=100*fuelcomp/sum(fuelcomp);
matcomp=[fuelcomp(1) fuelcomp(2) fuelcomp(1)+2*fuelcomp(2)+4*fuelcomp(3) fuelcomp(3)*UFrac];
saltzai=[30070 40090 90190];

MAT(1)=Mat('fuel'    ,1,-2.03,1.0050E+07,900,[saltzai UVec] ,matcomp);
MAT(2)=Mat('blanket' ,1,-4.44,1.9872E+07,900,[saltzai ThVec],[71 2 71+2*2+4*27 27*ThF]);
MAT(3)=Mat('graphite',2,-1.84,4.1807E+06,900,60000          ,1);
MAT(4)=Mat('feed'    ,0,0    ,1.0050E+07,900,UVec           ,10*MAT(1).atDens(MAT(1).find(UVec)));

volatiles=[2 10 18 36 44 45 46 47 49 54]; solubles=[1 7 8 34 35 41 42 43 52 53 38 39 56:1:65]; discard=[37 40 48 50 55];

REP(1)    =Rep('volatile1','fuel'   ,''       ,volatiles,[] ,1/60             ,'cont'); 
REP(end+1)=Rep('volatile2','blanket',''       ,volatiles,[] ,1/60             ,'cont'); 
REP(end+1)=Rep('solubles1','fuel'   ,''       ,solubles ,[] ,1/(110*24*3600)  ,'cont'); 
REP(end+1)=Rep('solubles2','blanket',''       ,solubles ,[] ,1/(110*24*3600)  ,'cont'); 
REP(end+1)=Rep('fuelproc' ,'blanket','fuel'   ,[92]     ,[] ,1/(1.1*24*3600)  ,'cont'); 
REP(end+1)=Rep('fuelproc' ,'blanket','feed'   ,[91]     ,[] ,1/(1.1*24*3600)  ,'cont');
REP(end+1)=Rep('discard1' ,'fuel'   ,''       ,discard  ,[] ,1/(5*365*24*3600),'cont'); 
REP(end+1)=Rep('discard2' ,'blanket',''       ,discard  ,[] ,1/(5*365*24*3600),'cont'); 
REP(end+1)=Rep('fuelproc' ,'feed'   ,'fuel'   ,922330   ,[] ,'keepAFPM'       ,'cont');
REP(end+1)=Rep('feed'     ,''       ,'blanket',ThVec    ,ThF,'keepAFPM'       ,'cont'); 
REP(end+1)=Rep('feed'     ,''       ,'blanket',ThVec    ,ThF,'keepAFPM'       ,'batch');
