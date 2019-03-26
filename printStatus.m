function [] = printStatus(MAT,SYS,prefix,time)
%PRINTSTATUS prints "burn" materials status to log
for j=SYS.IDX.MAT.burn
    fprintf('** %s ** %s : burn-up %.3f %%FIMA, %.3f %%mFP/(A+FP) at %s %02d cycle %02d time %G EFPD.\n',...
    prefix,MAT(j).name,MAT(j).FIMA,MAT(j).FPFrac,time,SYS.inCntr,SYS.ouCntr,SYS.nowTime(end));
end
return
end
