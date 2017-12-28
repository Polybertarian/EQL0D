function [] = saveFiles(MAT,OPT,SYS)

if(~OPT.PCC)
    suffix='';
else
    if(SYS.PCC.corrector)
        suffix='p';
    else
        suffix='c';
    end
end
if(SYS.ouCntr==1&~SYS.PCC.corrector)
    suffix='';
end
if(OPT.keepFiles)
    if(~SYS.stopOuter) % directory name
        SYS.dirName=['Cycle' num2str(SYS.ouCntr-1,'%03d') suffix];
    else
        SYS.dirName=['Cycle' num2str(SYS.ouCntr,'%03d') suffix];
    end
    mkdir(SYS.dirName);
    
    %%% move serpent output files to corresponding folder
    caseFiles={'_res.m','_arr0.m','_dep.m','_det0.m','.out'};
    for ext=caseFiles
        if(exist([SYS.Casename ext{1}],'file')==2)
            movefile([SYS.Casename ext{1}],SYS.dirName);
        end
    end
    for j=1:length(MAT)
        if(exist(['depmtx_' MAT(j).name '0.m'],'file')==2)
            movefile(['depmtx_' MAT(j).name '0.m'],SYS.dirName);
        end
        if(exist([MAT(j).name '.serp.bak'],'file')==2)
            movefile([MAT(j).name '.serp.bak'],SYS.dirName);
        end
    end

else
    if(exist([SYS.Casename '.out'],'file')==2)
        delete([SYS.Casename '.out']);
    end
end
    if(exist([SYS.Casename '.dep'],'file')==2)
        delete([SYS.Casename '.dep']);
    end
    if(exist('sss2.log','file')==2)
        delete('sss2.log');
    end
return
end
