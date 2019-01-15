function printStatus(MAT,RUN,prefix,time)
%PRINTSTATUS prints "burn" materials status to log

for j=1:numel(MAT) %SYS.IDX.MAT.burn
    fprintf('  ** %s **   %s : burn-up %.3f %%FIMA, %.3f %%mFP/(A+FP) at %s %02d cycle %02d time %G EFPD.\n',...
    prefix,MAT(j).name,MAT(j).FIMA,MAT(j).FPFrac,time,RUN.inCntr,RUN.ouCntr,RUN.nowTime(end));
end
return
end
