function nucName = ZAI2Name(nucList)
%ZAI2Name converts a vector of nuclide ZAIs to their name
if(~isempty(nucList))
        
listOfNuclides={'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',...
    'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',...
    'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo',...
    'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba',...
    'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',...
    'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',...
    'At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf',...
    'Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg'};

for i=1:length(nucList)
    if(nucList(i)<=111)
        if(nucList(i)~=-1)
            tmpStr=listOfNuclides{nucList(i)};
        else
            tmpStr='';
        end
    else
        if(nucList(i)~=666)
            stringZAI=strtrim(num2str(nucList(i)));
            metaIndicator=str2double(stringZAI(end));
            A=num2str(str2double(stringZAI(end-3:end-1)));
            if(strcmp(A,'0'))
                A='nat';
            end
            Z=str2double(stringZAI(1:end-4));
            tmpStr=[listOfNuclides{Z} '-' A];
            if(metaIndicator==1)
                tmpStr=[tmpStr 'm'];
            elseif(metaIndicator>1)
                tmpStr=[tmpStr 'm' num2str(metaIndicator)];
            end
        else
            tmpStr='';
        end
    end
    nucName{i}=tmpStr;   
end   
else
	nucName={};
end
nucName=nucName';
end