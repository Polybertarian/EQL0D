function coeffs = interpCoeffs(RUN)
    %INTERPCOEFFS gives predictor corrector interpolation coefficients
    if RUN.ouCntr==1&&~RUN.PCC.corrector||~RUN.PCC.active
        coeffs=[0 0 1]; %Constant extrapolation (no PCC or PCC@predictor@cycle 1)
    elseif RUN.ouCntr==1||RUN.PCC.active&&RUN.PCC.corrector
        coeffs=[0 1/2 1/2]; %linear interpolation (PCC@corrector@cycle 1)
    elseif RUN.ouCntr>1&&RUN.PCC.active&&~RUN.PCC.corrector
        coeffs=[0 -RUN.tStep(end)/(2*RUN.tStep(end-1)) 1+RUN.tStep(end)/(2*RUN.tStep(end-1))]; %linear extrapolation (PCC@predictor@cycle>1)
    else
        coeffs=[-RUN.tStep(end)^2/(6*RUN.tStep(end-1)*(RUN.tStep(end-1)*RUN.tStep(end))),...
        1/2+RUN.tStep(end)/(6*RUN.tStep(end-1)),...
        1/2-RUN.tStep/(6*(RUN.tStep(end)+RUN.tStep(end-1)))]; %quadratic interpolation (PCC@@cycle>1)
    end
    return
end
