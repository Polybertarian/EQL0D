load('endfb7')
OPT.nSteps=1;
%OPT.nCycles=35;
OPT.nCycles=1;
OPT.cycleLength=[1 3 3*[1 3 9 10 13 15 16 19 21 24 27 30 34 39 44 49 56 63 71 80 90]];
OPT.iterMode='steps';
OPT.reactControl=true;
OPT.writeMail=false;

ushare=0.112;
UVec=[942390 942400];
UFrac=[0.95 0.05];
OPT.REA.feedMat='feed';
OPT.REA.mode='addMass';
OPT.REA.upNuclides=UVec;
OPT.REA.downNuclides=[922330 922340];
OPT.REA.upFraction=UFrac;
OPT.REA.downFraction=[0.992 0.008];

%Thfrac=[0.0002 1-0.0002];
ThFrac=1;
%ThVec=[902300 902320];
ThVec=902320;

initfuel=[68.5 31.3 0.2]; fuelcomp=[initfuel(1:2) ushare]; fuelcomp=100*fuelcomp/sum(fuelcomp);
matcomp=[fuelcomp(1) fuelcomp(2) fuelcomp(1)+2*fuelcomp(2)+4*fuelcomp(3) fuelcomp(3)*UFrac];

MAT(1)=Mat('fuel'    ,1,-2.03,1.0050E+07,900,[30070 40090 90190 UVec],matcomp);
MAT(2)=Mat('blanket' ,1,-4.44,1.9872E+07,900,[30070 40090 90190 ThVec],[71.0 2.0 71.0+2*2.0+4*27.0 27*ThFrac]);
MAT(3)=Mat('graphite',2,-1.84,4.1807E+06,900,60000                    ,1);
MAT(4)=Mat('feed'    ,0,0    ,1.0050E+07,900,UVec                   ,10*MAT(1).atDens(MAT(1).find(UVec)));

volatiles=[2 10 18 36 44 45 46 47 49 54]; solubles=[1 7 8 34 35 41 42 43 52 53 38 39 56:1:65]; discard=[37 40 48 50 55];

REP(1).name='volatile';
REP(1).source='fuel';
REP(1).destination='void';
REP(1).elements=volatiles;
REP(1).share=1;
REP(1).rate=1/60;
REP(1).type='continuous';
REP(1).mode='remove';

REP(end+1).name='volatile2';
REP(end).source='blanket';
REP(end).destination='void';
REP(end).elements=volatiles;
REP(end).share=1;
REP(end).rate=1/60;
REP(end).type='continuous';
REP(end).mode='remove';

REP(end+1).name='fuelproc';
REP(end).source='fuel';
REP(end).destination='void';
REP(end).elements=solubles;
REP(end).share=1;
REP(end).rate=1/(110*24*3600);
REP(end).type='continuous';
REP(end).mode='remove';

REP(end+1).name='blanketproc';
REP(end).source='blanket';
REP(end).destination='void';
REP(end).elements=solubles;
REP(end).share=1;
REP(end).rate=1/(110*24*3600);
REP(end).type='continuous';
REP(end).mode='remove';

REP(end+1).name='storage2';
REP(end).source='blanket';
REP(end).destination='fuel';
REP(end).elements=92;
REP(end).share=1;
REP(end).rate=1/(1.1*24*3600);
REP(end).type='continuous';
REP(end).mode='remove';

REP(end+1).name='discard';
REP(end).source='fuel';
REP(end).destination='void';
REP(end).elements=discard;
REP(end).share=1;
REP(end).rate=1/(5*365*24*3600);
REP(end).type='continuous';
REP(end).mode='remove';

REP(end+1).name='discard2';
REP(end).source='blanket';
REP(end).destination='void';
REP(end).elements=discard;
REP(end).share=1;
REP(end).rate=1/(5*365*24*3600);
REP(end).type='continuous';
REP(end).mode='remove';

REP(end+1).name='refuel';
REP(end).source='void';
REP(end).destination='blanket';
REP(end).elements=ThVec;
REP(end).share=ThFrac;
REP(end).rate=1;
REP(end).type='continuous';
REP(end).mode='keepAFPM';
