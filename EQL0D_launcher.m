%%% Launcher for EQL0D
format long
warning('off','backtrace')
warning('off','MATLAB:declareGlobalBeforeUse')
warning('off','MATLAB:maxNumCompThreads:Deprecated')
warning('off','MATLAB:nearlySingularMatrix')

for i=1:length(Cases)
    if(exist(Cases{i},'file')==2)
        [PATH,FILE,EXT]=fileparts(Cases{i});
        if(isempty(PATH))
            SYS.Casename=[FILE EXT];
            EQL0D_main(SYS);
        else
            OLDDIR=cd(PATH);
            SYS.Casename=[FILE EXT];
            EQL0D_main(SYS);
            cd(OLDDIR);
        end 
    else
        warning(['Warning: file ' Cases{i} ' not found!'])
    end
end

