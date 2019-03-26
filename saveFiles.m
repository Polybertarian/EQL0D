function saveFiles(matNames,keepFiles,Casename,ouCntr,PCCActive)
    % SAVEFILES(matNames,keepFiles,SYS) moves Serpent outputs to be saved and deletes others

    outFiles={'depmtx_*.m','*_arr*.m','*_det*.m','*_dep.m','*.dep','*.out','*.bak'};

    if PCCActive % Delete all outputs if Predictor step
        keepFiles={};
    end

    files2delete=outFiles(~ismember(outFiles,keepFiles));
    for i=1:numel(files2delete)
        delete(files2delete{i})
    end
    caseFiles={'_arr','_det'};
    for ext=caseFiles
        if exist([Casename ext{1} '1.m'],'file')==2
            delete([Casename ext{1} '1.m']);
        end
    end

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
