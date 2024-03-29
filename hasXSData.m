function varargout=hasXSData(varargin)
%HASXSDATA Returns the ZAI given in varargin that have nuclear data libraries for neutron-induced reactions, discards the others

haveNuclearData=[
    10010	;	%	'H-1'
    10020	;	%	'H-2'
    10030	;	%	'H-3'
    20000	;	%	'He-nat'
    20030	;	%	'He-3'
    20040	;	%	'He-4'
    30060	;	%	'Li-6'
    30070	;	%	'Li-7'
    40070	;	%	'Be-7'
    40090	;	%	'Be-9'
    50100	;	%	'B-10'
    50110	;	%	'B-11'
    60000	;	%	'C-nat'
    60120	;	%	'C-12'
    70140	;	%	'N-14'
    70150	;	%	'N-15'
    80160	;	%	'O-16'
    80170	;	%	'O-17'
    90190	;	%	'F-19'
    100255	;	%	'Ne-25m5'
    110220	;	%	'Na-22'
    110230	;	%	'Na-23'
    120000	;	%	'Mg-nat'
    120240	;	%	'Mg-24'
    120250	;	%	'Mg-25'
    120260	;	%	'Mg-26'
    130270	;	%	'Al-27'
    140000	;	%	'Si-nat'
    140280	;	%	'Si-28'
    140290	;	%	'Si-29'
    140300	;	%	'Si-30'
    150310	;	%	'P-31'
    160000	;	%	'S-nat'
    160320	;	%	'S-32'
    160330	;	%	'S-33'
    160340	;	%	'S-34'
    160360	;	%	'S-36'
    170000	;	%	'Cl-nat'
    170350	;	%	'Cl-35'
    170370	;	%	'Cl-37'
    180360	;	%	'Ar-36'
    180380	;	%	'Ar-38'
    180400	;	%	'Ar-40'
    190000	;	%	'K-nat'
    190390	;	%	'K-39'
    190400	;	%	'K-40'
    190410	;	%	'K-41'
    200000	;	%	'Ca-nat'
    200400	;	%	'Ca-40'
    200420	;	%	'Ca-42'
    200430	;	%	'Ca-43'
    200440	;	%	'Ca-44'
    200460	;	%	'Ca-46'
    200480	;	%	'Ca-48'
    210450	;	%	'Sc-45'
    220000	;	%	'Ti-nat'
    220460	;	%	'Ti-46'
    220470	;	%	'Ti-47'
    220480	;	%	'Ti-48'
    220490	;	%	'Ti-49'
    220500	;	%	'Ti-50'
    230000	;	%	'V-nat'
    230510	;	%	'V-51'
    240000	;	%	'Cr-nat'
    240500	;	%	'Cr-50'
    240520	;	%	'Cr-52'
    240530	;	%	'Cr-53'
    240540	;	%	'Cr-54'
    250550	;	%	'Mn-55'
    260000	;	%	'Fe-nat'
    260540	;	%	'Fe-54'
    260560	;	%	'Fe-56'
    260570	;	%	'Fe-57'
    260580	;	%	'Fe-58'
    270580	;	%	'Co-58'
    270590	;	%	'Co-59'
    271580	;	%	'Co-158'
    271581	;	%	'Co-158m'
    280000	;	%	'Ni-nat'
    280580	;	%	'Ni-58'
    280590	;	%	'Ni-59'
    280600	;	%	'Ni-60'
    280610	;	%	'Ni-61'
    280620	;	%	'Ni-62'
    280640	;	%	'Ni-64'
    290000	;	%	'Cu-nat'
    290630	;	%	'Cu-63'
    290650	;	%	'Cu-65'
    300000	;	%	'Zn-nat'
    300640	;	%	'Zn-64'
    310000	;	%	'Ga-nat'
    310690	;	%	'Ga-69'
    310710	;	%	'Ga-71'
    320000	;	%	'Ge-nat'
    320700	;	%	'Ge-70'
    320720	;	%	'Ge-72'
    320730	;	%	'Ge-73'
    320740	;	%	'Ge-74'
    320760	;	%	'Ge-76'
    330740	;	%	'As-74'
    330750	;	%	'As-75'
    340740	;	%	'Se-74'
    340760	;	%	'Se-76'
    340770	;	%	'Se-77'
    340780	;	%	'Se-78'
    340790	;	%	'Se-79'
    340800	;	%	'Se-80'
    340820	;	%	'Se-82'
    350790	;	%	'Br-79'
    350810	;	%	'Br-81'
    360780	;	%	'Kr-78'
    360800	;	%	'Kr-80'
    360820	;	%	'Kr-82'
    360830	;	%	'Kr-83'
    360840	;	%	'Kr-84'
    360850	;	%	'Kr-85'
    360860	;	%	'Kr-86'
    370850	;	%	'Rb-85'
    370860	;	%	'Rb-86'
    370870	;	%	'Rb-87'
    380840	;	%	'Sr-84'
    380860	;	%	'Sr-86'
    380870	;	%	'Sr-87'
    380880	;	%	'Sr-88'
    380890	;	%	'Sr-89'
    380900	;	%	'Sr-90'
    390890	;	%	'Y-89'
    390900	;	%	'Y-90'
    390910	;	%	'Y-91'
    400000	;	%	'Zr-nat'
    400900	;	%	'Zr-90'
    400910	;	%	'Zr-91'
    400920	;	%	'Zr-92'
    400930	;	%	'Zr-93'
    400940	;	%	'Zr-94'
    400950	;	%	'Zr-95'
    400960	;	%	'Zr-96'
    410930	;	%	'Nb-93'
    410940	;	%	'Nb-94'
    410950	;	%	'Nb-95'
    420000	;	%	'Mo-nat'
    420920	;	%	'Mo-92'
    420940	;	%	'Mo-94'
    420950	;	%	'Mo-95'
    420960	;	%	'Mo-96'
    420970	;	%	'Mo-97'
    420980	;	%	'Mo-98'
    420990	;	%	'Mo-99'
    421000	;	%	'Mo-100'
    430990	;	%	'Tc-99'
    440960	;	%	'Ru-96'
    440980	;	%	'Ru-98'
    440990	;	%	'Ru-99'
    441000	;	%	'Ru-100'
    441010	;	%	'Ru-101'
    441020	;	%	'Ru-102'
    441030	;	%	'Ru-103'
    441040	;	%	'Ru-104'
    441050	;	%	'Ru-105'
    441060	;	%	'Ru-106'
    451030	;	%	'Rh-103'
    451050	;	%	'Rh-105'
    461020	;	%	'Pd-102'
    461040	;	%	'Pd-104'
    461050	;	%	'Pd-105'
    461060	;	%	'Pd-106'
    461070	;	%	'Pd-107'
    461080	;	%	'Pd-108'
    461100	;	%	'Pd-110'
    470000	;	%	'Ag-nat'
    471070	;	%	'Ag-107'
    471090	;	%	'Ag-109'
    471101	;	%	'Ag-110m'
    471110	;	%	'Ag-111'
    480000	;	%	'Cd-nat'
    481060	;	%	'Cd-106'
    481080	;	%	'Cd-108'
    481100	;	%	'Cd-110'
    481110	;	%	'Cd-111'
    481120	;	%	'Cd-112'
    481130	;	%	'Cd-113'
    481140	;	%	'Cd-114'
    481151	;	%	'Cd-115m'
    481160	;	%	'Cd-116'
    490000	;	%	'In-nat'
    491130	;	%	'In-113'
    491150	;	%	'In-115'
    500000	;	%	'Sn-nat'
    501120	;	%	'Sn-112'
    501130	;	%	'Sn-113'
    501140	;	%	'Sn-114'
    501150	;	%	'Sn-115'
    501160	;	%	'Sn-116'
    501170	;	%	'Sn-117'
    501180	;	%	'Sn-118'
    501190	;	%	'Sn-119'
    501200	;	%	'Sn-120'
    501220	;	%	'Sn-122'
    501230	;	%	'Sn-123'
    501240	;	%	'Sn-124'
    501250	;	%	'Sn-125'
    501260	;	%	'Sn-126'
    510000	;	%	'Sb-nat'
    511210	;	%	'Sb-121'
    511230	;	%	'Sb-123'
    511240	;	%	'Sb-124'
    511250	;	%	'Sb-125'
    511260	;	%	'Sb-126'
    521200	;	%	'Te-120'
    521220	;	%	'Te-122'
    521230	;	%	'Te-123'
    521240	;	%	'Te-124'
    521250	;	%	'Te-125'
    521260	;	%	'Te-126'
    521271	;	%	'Te-127m'
    521280	;	%	'Te-128'
    521291	;	%	'Te-129m'
    521300	;	%	'Te-130'
    521320	;	%	'Te-132'
    531270	;	%	'I-127'
    531290	;	%	'I-129'
    531300	;	%	'I-130'
    531310	;	%	'I-131'
    531350	;	%	'I-135'
    541230	;	%	'Xe-123'
    541240	;	%	'Xe-124'
    541260	;	%	'Xe-126'
    541280	;	%	'Xe-128'
    541290	;	%	'Xe-129'
    541300	;	%	'Xe-130'
    541310	;	%	'Xe-131'
    541320	;	%	'Xe-132'
    541330	;	%	'Xe-133'
    541340	;	%	'Xe-134'
    541350	;	%	'Xe-135'
    541360	;	%	'Xe-136'
    551330	;	%	'Cs-133'
    551340	;	%	'Cs-134'
    551350	;	%	'Cs-135'
    551360	;	%	'Cs-136'
    551370	;	%	'Cs-137'
    561300	;	%	'Ba-130'
    561320	;	%	'Ba-132'
    561330	;	%	'Ba-133'
    561340	;	%	'Ba-134'
    561350	;	%	'Ba-135'
    561360	;	%	'Ba-136'
    561370	;	%	'Ba-137'
    561380	;	%	'Ba-138'
    561400	;	%	'Ba-140'
    571380	;	%	'La-138'
    571390	;	%	'La-139'
    571400	;	%	'La-140'
    581360	;	%	'Ce-136'
    581380	;	%	'Ce-138'
    581390	;	%	'Ce-139'
    581400	;	%	'Ce-140'
    581410	;	%	'Ce-141'
    581420	;	%	'Ce-142'
    581430	;	%	'Ce-143'
    581440	;	%	'Ce-144'
    591410	;	%	'Pr-141'
    591420	;	%	'Pr-142'
    591430	;	%	'Pr-143'
    601420	;	%	'Nd-142'
    601430	;	%	'Nd-143'
    601440	;	%	'Nd-144'
    601450	;	%	'Nd-145'
    601460	;	%	'Nd-146'
    601470	;	%	'Nd-147'
    601480	;	%	'Nd-148'
    601500	;	%	'Nd-150'
    611470	;	%	'Pm-147'
    611480	;	%	'Pm-148'
    611481  ;   %   'Pm-148m'
    611490	;	%	'Pm-149'
    611510	;	%	'Pm-151'
    620000	;	%	'Sm-nat'
    621440	;	%	'Sm-144'
    621470	;	%	'Sm-147'
    621480	;	%	'Sm-148'
    621490	;	%	'Sm-149'
    621500	;	%	'Sm-150'
    621510	;	%	'Sm-151'
    621520	;	%	'Sm-152'
    621530	;	%	'Sm-153'
    621540	;	%	'Sm-154'
    630000	;	%	'Eu-nat'
    631510	;	%	'Eu-151'
    631520	;	%	'Eu-152'
    631530	;	%	'Eu-153'
    631540	;	%	'Eu-154'
    631550	;	%	'Eu-155'
    631560	;	%	'Eu-156'
    631570	;	%	'Eu-157'
    640000	;	%	'Gd-nat'
    641520	;	%	'Gd-152'
    641530	;	%	'Gd-153'
    641540	;	%	'Gd-154'
    641550	;	%	'Gd-155'
    641560	;	%	'Gd-156'
    641570	;	%	'Gd-157'
    641580	;	%	'Gd-158'
    641600	;	%	'Gd-160'
    651590	;	%	'Tb-159'
    651600	;	%	'Tb-160'
    661560	;	%	'Dy-156'
    661580	;	%	'Dy-158'
    661600	;	%	'Dy-160'
    661610	;	%	'Dy-161'
    661620	;	%	'Dy-162'
    661630	;	%	'Dy-163'
    661640	;	%	'Dy-164'
    671650	;	%	'Ho-165'
    671661	;	%	'Ho-166m'
    681620	;	%	'Er-162'
    681640	;	%	'Er-164'
    681660	;	%	'Er-166'
    681670	;	%	'Er-167'
    681680	;	%	'Er-168'
    681700	;	%	'Er-170'
    710000	;	%	'Lu-nat'
    711750	;	%	'Lu-175'
    711760	;	%	'Lu-176'
    720000	;	%	'Hf-nat'
    721740	;	%	'Hf-174'
    721760	;	%	'Hf-176'
    721770	;	%	'Hf-177'
    721780	;	%	'Hf-178'
    721790	;	%	'Hf-179'
    721800	;	%	'Hf-180'
    731810	;	%	'Ta-181'
    731820	;	%	'Ta-182'
    740000	;	%	'W-nat'
    741820	;	%	'W-182'
    741830	;	%	'W-183'
    741840	;	%	'W-184'
    741860	;	%	'W-186'
    750000	;	%	'Re-nat'
    751850	;	%	'Re-185'
    751870	;	%	'Re-187'
    760000	;	%	'Os-nat'
    770000	;	%	'Ir-nat'
    771910	;	%	'Ir-191'
    771930	;	%	'Ir-193'
    780000	;	%	'Pt-nat'
    791970	;	%	'Au-197'
    800000	;	%	'Hg-nat'
    801960	;	%	'Hg-196'
    801980	;	%	'Hg-198'
    801990	;	%	'Hg-199'
    802000	;	%	'Hg-200'
    802010	;	%	'Hg-201'
    802020	;	%	'Hg-202'
    802040	;	%	'Hg-204'
    810000	;	%	'Tl-nat'
    820000	;	%	'Pb-nat'
    822040	;	%	'Pb-204'
    822060	;	%	'Pb-206'
    822070	;	%	'Pb-207'
    822080	;	%	'Pb-208'
    832090	;	%	'Bi-209'
    882230	;	%	'Ra-223'
    882240	;	%	'Ra-224'
    882250	;	%	'Ra-225'
    882260	;	%	'Ra-226'
    892250	;	%	'Ac-225'
    892260	;	%	'Ac-226'
    892270	;	%	'Ac-227'
    902270	;	%	'Th-227'
    902280	;	%	'Th-228'
    902290	;	%	'Th-229'
    902300	;	%	'Th-230'
    902320	;	%	'Th-232'
    902330	;	%	'Th-233'
    902340	;	%	'Th-234'
    912310	;	%	'Pa-231'
    912320	;	%	'Pa-232'
    912330	;	%	'Pa-233'
    922320	;	%	'U-232'
    922330	;	%	'U-233'
    922340	;	%	'U-234'
    922350	;	%	'U-235'
    922360	;	%	'U-236'
    922370	;	%	'U-237'
    922380	;	%	'U-238'
    922390	;	%	'U-239'
    922400	;	%	'U-240'
    922410	;	%	'U-241'
    932350	;	%	'Np-235'
    932360	;	%	'Np-236'
    932370	;	%	'Np-237'
    932380	;	%	'Np-238'
    932390	;	%	'Np-239'
    942360	;	%	'Pu-236'
    942370	;	%	'Pu-237'
    942380	;	%	'Pu-238'
    942390	;	%	'Pu-239'
    942400	;	%	'Pu-240'
    942410	;	%	'Pu-241'
    942420	;	%	'Pu-242'
    942430	;	%	'Pu-243'
    942440	;	%	'Pu-244'
    942460	;	%	'Pu-246'
    952410	;	%	'Am-241'
    952420	;	%	'Am-242'
    952421	;	%	'Am-242m'
    952430	;	%	'Am-243'
    952440	;	%	'Am-244'
    952441  ;   %   'Am-244m'
    962400	;	%	'Cm-240'
    962410	;	%	'Cm-241'
    962420	;	%	'Cm-242'
    962430	;	%	'Cm-243'
    962440	;	%	'Cm-244'
    962450	;	%	'Cm-245'
    962460	;	%	'Cm-246'
    962470	;	%	'Cm-247'
    962480	;	%	'Cm-248'
    962490	;	%	'Cm-249'
    962500	;	%	'Cm-250'
    972470	;	%	'Bk-247'
    972490	;	%	'Bk-249'
    972500	;	%	'Bk-250'
    982490	;	%	'Cf-249'
    982500	;	%	'Cf-250'
    982510	;	%	'Cf-251'
    982520	;	%	'Cf-252'
    982530	;	%	'Cf-253'
    982540	;	%	'Cf-254'
    992530	;	%	'Es-253'
    992540	;	%	'Es-254'
    992550		%	'Es-255'
    ]';

varargout = {ismember([varargin{:}],haveNuclearData)};
return
end

