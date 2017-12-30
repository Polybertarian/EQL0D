function MAT = redoxControl(MAT,OPT,SYS)
%REDOXCONTROL Controls concentration of F/Cl ions

for i=SYS.IDX.redoxMat
    default=MAT(i).halideExcess*MAT(i).volume;
    switch OPT.REDOX.replaceMode
        case 'replaceMass'
            dN=default/(1+MAT(i).oxState(SYS.IDX.redoxNuc{i})*MAT(i).avMass(SYS.IDX.redoxHalide{i})/...
                MAT(i).atomicMass(SYS.IDX.redoxNuc{i}));
            MAT(i).N(SYS.IDX.redoxNuc{i},end+1)=MAT(i).N(SYS.IDX.redoxNuc{i},end)-dN*...
                MAT(i).avMass(SYS.IDX.redoxHalide{i})/MAT(i).atomicMass(SYS.IDX.redoxNuc{i});
        case 'replaceAtom'
            dN=default/(1+MAT(i).oxState(SYS.IDX.redoxNuc{i}));
            MAT(i).N(SYS.IDX.redoxNuc{i},end+1)=MAT(i).N(SYS.IDX.redoxNuc{i},end)-dN;
        otherwise
            dN=default;
    end
    MAT(i).N(SYS.IDX.redoxHalide{i},end+1)=MAT(i).N(SYS.IDX.redoxHalide{i},end)+dN*MAT(i).aFrac(SYS.IDX.redoxHalide{i});
    name=MAT.nuclideName(SYS.IDX.redoxHalide{i});
    fprintf(SYS.FID.log,'%s\n',['** REDOX ** Excess of ' ...
        num2str(-default*1E24,'%E') ' ' [name{:}]  ' atoms corrected.']);   
end
return
end

