function coeffs = interpCoeffs(SYS)
%INTERPCOEFFS gives predictor corrector interpolation coefficients 
if(SYS.ouCntr==1&&~SYS.PCC.corrector||~SYS.PCC.active)
    coeffs=[0 0 1]; %Constant extrapolation (no PCC or PCC@predictor@cycle 1)
elseif(SYS.ouCntr==1||SYS.PCC.active&&SYS.PCC.corrector)
    coeffs=[0 1/2 1/2]; %linear interpolation (PCC@corrector@cycle 1)
elseif(SYS.ouCntr>1&&SYS.PCC.active&&~SYS.PCC.corrector)
    coeffs=[0 -SYS.tStep(end)/(2*SYS.tStep(end-1)) 1+SYS.tStep(end)/(2*SYS.tStep(end-1))]; %linear extrapolation (PCC@predictor@cycle>1)
else
    coeffs=[-SYS.tStep(end)^2/(6*SYS.tStep(end-1)*(SYS.tStep(end-1)*SYS.tStep(end))),...
        1/2+SYS.tStep(end)/(6*SYS.tStep(end-1)),...
        1/2-SYS.tStep/(6*(SYS.tStep(end)+SYS.tStep(end-1)))]; %quadratic interpolation (PCC@@cycle>1)
end
return
end