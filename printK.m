function [ ] = printK(SYS, step, prefix, source)
%PRINTK=(SYS,step,prefix,source) Prints K-eff values to keff.txt

global FID

switch source
  case 'EQL0D'
    fprintf(FID.keff,'%-7s %-3s %5.3d %4.3d %3s %-12.6G %-9f %-9f\n',...
      'EQL0D',prefix,SYS.RUN.ouCntr,SYS.RUN.inCntr,step,SYS.RUN.nowTime(end),SYS.KINF.EQL0D(end),SYS.KEFF.EQL0D(end));
  case 'Serpent'
    fprintf(FID.keff,'%-7s %-3s %5.3d %4.3d %3s %-12.6G %-9f %-9f\n',...
      'Serpent',prefix,SYS.RUN.ouCntr,SYS.RUN.inCntr,step,SYS.RUN.nowTime(end),SYS.KINF.Serpent(end),SYS.KEFF.Serpent(end));
end

return
end