function [ ] = printK(SYS, step, prefix, source)
%PRINTK=(SYS,step,prefix,source) Prints K-eff values to keff.txt

switch source
  case 'EQL0D'
    fprintf(SYS.FID.keff,'%-7s %-3s %5.3d %4.3d %3s %-12.6G %-9f %-9f\n',...
      'EQL0D',prefix,SYS.ouCntr,SYS.inCntr,step,SYS.nowTime(end),SYS.KINF.EQL0D(end),SYS.KEFF.EQL0D(end));
  case 'Serpent'
    fprintf(SYS.FID.keff,'%-7s %-3s %5.3d %4.3d %3s %-12.6G %-9f %-9f\n',...
      'Serpent',prefix,SYS.ouCntr,SYS.inCntr,step,SYS.nowTime(end),SYS.KINF.Serpent(end),SYS.KEFF.Serpent(end));
end

return
end