function [] = modifyInput(inputName,paramName,paramValue)
%MODIFYINPUT changes the parameter paramName to the new value paramValue in
%Serpent input file inputName

% Adapt to different versions of 'grep'
if(~ismac)
    pos=1;
else
    pos=2;
end

if(exist(inputName,'file')==2)
    [~,output]=unix(['grep -nr ''' paramName ' '' ' inputName]);
    output=textscan(output,'%s','Delimiter',{' ',':'});
    line=output{1}(pos);
    switch paramName
        case 'dep'
            currentValue=output{1}(pos+3);
            output=strcat('sed -i.bak ''',line,'s/',currentValue,'/', num2str(paramValue),'/''',[' ' inputName]);
            [status,~]=unix(output{1});
        case 'depmtx'
            currentValue=output{1}(pos+3);
            output=strcat('sed -i.bak ''',line,'s/',currentValue,'/', num2str(paramValue),'/g''',[' ' inputName]);
            [status,~]=unix(output{1});
        otherwise
            error(['Error: parameter to change ''' paramName ''' not recognized!'])
    end
    if(status~=0)
        error('Error while writing in Serpent file!');
    end
end

end