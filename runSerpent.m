function runSerpent(serpentPath,Casename,nCores)
%RUNSERPENT Calls serpentPath with input fileName and nCores

if exist([Casename '.seed'],'file')==2
    [status,errmsg]=unix([serpentPath ' -omp ' num2str(nCores) ' -replay ' Casename ' > sss2.log']);
else
    [status,errmsg]=unix([serpentPath ' -omp ' num2str(nCores) ' ' Casename ' > sss2.log']);
end

if status~=0
    error('Serpent:CrashDuringRun',['Error: SSS simulation aborted \n' errmsg]);
else
    fprintf('%s\n','*** SERPENT *** Serpent run finished!');
end
return
end
