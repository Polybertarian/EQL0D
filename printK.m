function [ ] = printK(SYS, step, prefix, source)
%PRINTK=(SYS,step,prefix,source) Prints K-eff values to keff.txt

fprintf(SYS.FID.keff,'%-7s %-3s %5.3d %4.3d %3s %-12.6G %-9f %-9f\n',...
    source,prefix,SYS.ouCntr,SYS.inCntr,step,SYS.nowTime(end),SYS.kinf(end),SYS.keff(end));

return
end