function EQL0Dlauncher(Cases,nCores,verbose,printNquit,restart,reset,debug)
    % Launcher for EQL0D, handles batch case launching
    SYS=struct('nCores',nCores,'verboseMode',logical(verbose),'printAndQuit',logical(printNquit),...
    'restartCalc',logical(restart),'resetCounters',logical(reset),'debugMode',logical(debug));
    format long
    warning('off','backtrace')
    warning('off','MATLAB:declareGlobalBeforeUse')
    warning('off','MATLAB:maxNumCompThreads:Deprecated')
    warning('off','MATLAB:nearlySingularMatrix')

    for i=1:length(Cases)
        if exist(Cases{i},'file')==2
            [PATH,FILE,EXT]=fileparts(Cases{i});
            if strcmp(EXT,'.m')
                EXT='';
            end
            SYS.Casename=[FILE EXT];
            if ~isempty(PATH)
                OLDDIR=cd(PATH);
            end
            if ~restart&~reset
                if numel([dir('depmtx_*.m');dir('*_res*.m');dir('*_arr*.m');dir('*_det*.m');dir('*.mat');dir('*.txt');dir('*.log');])>0&~debug
                    j=1;
                    while numel(dir(['run' num2str(j,'%.2d')]))>0
                        j=j+1;
                    end
                    dirName=['run' num2str(j,'%.2d')];
                    mkdir(dirName);
                    files2Move=[dir('depmtx_*.m');dir('*_res*.m');dir('*_arr*.m');dir('*_det*.m');dir('*.mat');dir('*.txt');dir([SYS.Casename '.log'])];
                    files2Move={files2Move.name};
                    for k=1:numel(files2Move)
                        movefile(files2Move{k},dirName);
                    end
                end
            end
            main(SYS);
            if ~isempty(PATH)
                cd(OLDDIR);
            end
        else
            warning(['Warning: file ' Cases{i} ' not found! Skipping...'])
        end
    end
end
