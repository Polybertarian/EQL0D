function EQL0Dlauncher(Cases,nCores,verbose,printNquit,restart,reset,debug)
%%% Launcher for EQL0D
SYS=struct('nCores',nCores,'verboseMode',logical(verbose),'printAndQuit',logical(printNquit),...
    'restartCalc',logical(restart),'resetCounters',logical(reset),'debugMode',logical(debug));
format long 
warning('off','backtrace')
warning('off','MATLAB:declareGlobalBeforeUse')
warning('off','MATLAB:maxNumCompThreads:Deprecated')
warning('off','MATLAB:nearlySingularMatrix')

for i=1:length(Cases)
    if(exist(Cases{i},'file')==2)
        [PATH,FILE,EXT]=fileparts(Cases{i});
        if(strcmp(EXT,'.m'))
          EXT='';
        end
        if(isempty(PATH))
            SYS.Casename=[FILE EXT];
            EQL0D(SYS);
        else
            OLDDIR=cd(PATH);
            SYS.Casename=[FILE EXT];
            EQL0D(SYS);
            cd(OLDDIR);
        end 
    else
        warning(['Warning: file ' Cases{i} ' not found! Skipping...'])
    end
end

end