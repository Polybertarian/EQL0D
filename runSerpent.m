function [] = runSerpent(OPT,SYS)
%RUNSERPENT Calls OPT.serpentPath with input fileName and SYS.nCores
if(SYS.verboseMode)
    fprintf(SYS.FID.log,'%s\n','*** SERPENT *** Running Serpent...');
end
if(exist([SYS.Casename '.seed'],'file')==2)
    [status,~]=unix([OPT.serpentPath ' -omp ' num2str(SYS.nCores) ' -replay ' SYS.Casename ' > sss2.log']);
else
	[status,~]=unix([OPT.serpentPath ' -omp ' num2str(SYS.nCores) ' ' SYS.Casename ' > sss2.log']);
end

if(status~=0)
	error('Error: SSS simulation aborted');
elseif(SYS.verboseMode)
    fprintf(SYS.FID.log,'%s\n','*** SERPENT *** Serpent run finished!');
end
return
end

