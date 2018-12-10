function openLogs(restartCalc,REA)
%OPENLOGS Opens log files for EQL0D

global FID

if restartCalc
    permission='at';
else
    permission='wt';
end

[FID.keff,errmsg]=fopen('keff.txt',permission);
if ~restartCalc
    fprintf(FID.keff,'%-7s %-3s %-5s %-4s %-3s %-12s %-9s %-9s\n','Source','PCC','Cycle','Step','Rep','Time','k-inf','k-eff');
end

if ~isempty(REA)
    switch REA.mode
        case {'replace','addMass'}
            [FID.react,errmsg]=fopen('reactivity.txt',permission);
            if ~restartCalc
                name{1}=REA.targetMat;
                if ~isempty(REA.feedMat)
                    name{2}=REA.feedMat;
                else
                    name{2}=[];
                end
                fprintf(FID.react,'%-22s','Time');
                elementsTarget=unique([REA.replNuclides;REA.upNuclides;REA.downNuclides]);
                fprintf(FID.react,['%-' num2str(13*numel(elementsTarget)) 's'],name{1});
                if ~isempty(REA.feedMat)
                    fprintf(FID.react,['%-' num2str(13*numel(elementsTarget)) 's'],name{2});
                end
                fprintf(FID.react,'\n');
                fprintf(FID.react,'%-7s%-6s%-9s','Cycle','Step','EFPD');
                names{1}=ZAI2Name(elementsTarget);
                fprintf(FID.react,repmat('%-13s',1,numel(names{1})),names{1}{:});
                if ~isempty(REA.feedMat)
                    fprintf(FID.react,repmat('%-13s',1,numel(names{1})),names{1}{:});
                end
                fprintf(FID.react,'\n');
            end
        case 'addVolume'
            [FID.volume,errmsg]=fopen('volume.txt',permission);
            if ~restartCalc
                fprintf(FID.volume,'%-7s%-6s%-9s%-12s\n','Cycle','Step','EFPD','Added Volume [cm^3]');
            end
    end
end

end

