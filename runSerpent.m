function [] = runSerpent(serpentPath,Casename,nCores,logFID)
%RUNSERPENT Calls OPT.serpentPath with input fileName and nCores

if(exist([Casename '.seed'],'file')==2)
    [status,~]=unix([serpentPath ' -omp ' num2str(nCores) ' -replay ' Casename ' > sss2.log']);
else
    [status,~]=unix([serpentPath ' -omp ' num2str(nCores) ' ' Casename ' > sss2.log']);
end

if(status~=0)
    error('Error: SSS simulation aborted');
else
    fprintf(logFID,'%s\n','*** SERPENT *** Serpent run finished!');
end
return
end

