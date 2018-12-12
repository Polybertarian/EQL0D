function saveFiles(matNames,keepFiles,Casename,ouCntr,PCCActive)
% SAVEFILES(matNames,keepFiles,SYS) moves Serpent outputs to be saved and deletes others
caseFiles={'_arr','_det'};
for ext=caseFiles
    if exist([Casename ext{1} '1.m'],'file')==2
       delete([Casename ext{1} '1.m']);
    end
end
delete('*_dep.m'); delete('*.dep'); delete('sss2.log'); delete('*.out'); delete('*.bak');
if ~keepFiles||PCCActive
    delete('depmtx_*.m'); delete('*_det*.m'); delete('*_arr*.m');
else
    cycleExt=num2str(ouCntr-1,'%.2d');
    if ouCntr>1
        if exist([Casename '_res' '.m'],'file')==2
            movefile([Casename '_res' '.m'],[Casename '_res' cycleExt '.m']);
        end
        for ext=caseFiles
            if exist([Casename ext{1} '0.m'],'file')==2
                movefile([Casename ext{1} '0.m'],[Casename ext{1} cycleExt '.m']);
            end
        end
        for j=1:length(matNames)
            if exist(['depmtx_' matNames{j} '0.m'],'file')==2
                movefile(['depmtx_' matNames{j} '0.m'],['depmtx_' matNames{j} cycleExt '.m']);
            end
        end
    end
end
return
end
