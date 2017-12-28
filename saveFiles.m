function [] = saveFiles(matNames,keepFiles,SYS)
% SAVEFILES(matNames,keepFiles,SYS) moves Serpent outputs to be saved to corresponding folder and deletes
% others

if(keepFiles)
    if(~SYS.PCC.active)
        dirName=['Cycle' num2str(SYS.ouCntr-1,'%03d')];
    else
        if(SYS.PCC.corrector)
            dirName=['Cycle' num2str(SYS.ouCntr,'%03d') 'p'];
        else
            if(SYS.ouCntr==1)
                dirName=['Cycle' num2str(SYS.ouCntr-1,'%03d')];
            else
                dirName=['Cycle' num2str(SYS.ouCntr-1,'%03d') 'c'];
            end
        end
    end
    mkdir(dirName);
    
    %%% move serpent output files to corresponding folder
    caseFiles={'_res.m','_arr0.m','_dep.m','_det0.m','.out'};
    for ext=caseFiles
        if(exist([SYS.Casename ext{1}],'file')==2)
            movefile([SYS.Casename ext{1}],dirName);
        end
    end
    for j=1:length(matNames)
        if(exist(['depmtx_' matNames{j} '0.m'],'file')==2)
            movefile(['depmtx_' matNames{j} '0.m'],dirName);
        end
        if(exist([matNames{j} '.serp.bak'],'file')==2)
            movefile([matNames{j} '.serp.bak'],dirName);
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
