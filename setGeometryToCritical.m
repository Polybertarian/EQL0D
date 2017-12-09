function status = setGeometryToCritical(OPT)
%SETGEOMETRYTOCRITICAL Changes the geometry to the critical value for a cylinder.
SYS.Casename

% Adapt to different versions of 'grep'
if(~ismac)
    pos=1;
else
    pos=2;
end

run([Casename '_res.m']);

if(INF_KINF(end,1)>1)
    %%% Calculate critical reactor radius
    INF_MIGRATION_AREA=INF_DIFFCOEF(end,[1 3])./[INF_REMXS(end,1) INF_ABS(end,3)];
    INF_MIGRATION_AREA(isnan(INF_MIGRATION_AREA))=[];
    INF_MIGRATION_AREA=sum(INF_MIGRATION_AREA);
    INF_CRITICAL_RADIUS=2.405*sqrt(3/2)*sqrt(INF_MIGRATION_AREA/(INF_KINF(end,1)-1));

    %%% Change in geometry.serp
    [~,output]=unix('grep -nr "surf core cyl" geometry.serp');
    output=textscan(output,'%s','Delimiter',{' ',':'})
    line=output{1}(pos)
    currentValue=output{1}(pos+6);
    output=strcat('sed -i.bak ''',line,'s/',currentValue,'/', num2str(INF_CRITICAL_RADIUS),'/g''',' geometry.serp');
    [~,~]=unix(output{1});
else
    warning('k-inf < 1: Impossible to make finite geometry critical! Skipping...')
end

end
