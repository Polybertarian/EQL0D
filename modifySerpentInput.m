function modifySerpentInput(inputName,paramName,paramValue)
%MODIFYINPUT changes the parameter paramName to the new value paramValue in
%Serpent input file inputName

% Adapt to different versions of grep
pos=1+ismac;

if(exist(inputName,'file')==2)
    [~,output]=unix(['grep -nr ''' paramName ' '' ' inputName]);
    output=textscan(output,'%s','Delimiter',{' ',':'});
    line=output{1}(pos);
    switch paramName
        case 'dep'
            currentValue=output{1}(pos+3);
            output=strcat('sed -i.bak ''',line,'s/',currentValue,'/', num2str(paramValue),'/''',[' ' inputName]);
            [status,errmsg]=unix(output{1});
        case 'det'
            isAbsent=[]; % Check presence of necessary input parameters in Serpent file
            %[~,isAbsent{end+1}]=unix(['grep -c "set depmtx 1" ' inputName]);
            [~,isAbsent{end+1}]=unix(['grep -c "set arr 1" ' inputName]);
            [~,isAbsent{end+1}]=unix(['grep -c "det intFlux" ' inputName]);
            [~,isAbsent{end+1}]=unix(['grep -c "det intFiss" ' inputName]);
            [~,isAbsent{end+1}]=unix(['grep -c "det intProd" ' inputName]);
            [~,isAbsent{end+1}]=unix(['grep -c "det intCapt" ' inputName]);
            isAbsent=~logical(str2double(isAbsent));
            
            if(any(isAbsent))
                [fID,~]=fopen(inputName,'a');
                % if(isAbsent(1))
                % fprintf(fID,'\n%s\n','set depmtx 1 1');
                % end
                if(isAbsent(1))
                    fprintf(fID,'\n%s\n','set arr 1 0');
                end
                if(isAbsent(2))
                    fprintf(fID,'\n%s','det intFlux ');
                    for i=SYS.IDX.MAT.inFlux
                        fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
                    end
                end
                if(isAbsent(3))
                    fprintf(fID,'\n%s','det intFiss dr -6 void ');
                    for i=SYS.IDX.MAT.inFlux
                        fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
                    end
                end
                if(isAbsent(4))
                    fprintf(fID,'\n%s','det intProd dr -7 void ');
                    for i=SYS.IDX.MAT.inFlux
                        fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
                    end
                end
                if(isAbsent(5))
                    fprintf(fID,'\n%s','det intCapt dr -2 void ');
                    for i=SYS.IDX.MAT.inFlux
                        fprintf(fID,'%s',['dm ' MAT(i).name ' ']);
                    end
                end
                fclose(fID);
            end
        otherwise
            error(['Error: parameter to change ''' paramName ''' not recognized!'])
    end
    if(status~=0)
        error(['Error: while writing in input: ' errmsg]);
    end
else
    error('Error: Input file not found!')
end

end
