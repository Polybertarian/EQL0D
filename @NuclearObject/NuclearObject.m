classdef NuclearObject
  %NUCLEAROBJECT Object that has nuclear properties (half-life, atomic
  %mass, etc...
  
  properties
    atomicMass
    decayEnergy
    halfLife
    ZAI
    oxState
    ingToxicity
    inhToxicity
    hasNucData
    nuclideName
    defDecMtx
    burnIdx
  end
  
  methods
    function obj = NuclearObject()
      global DAT
      if exist('DAT','var')
        obj.burnIdx=isProduced(DAT.libraryName,DAT.ZAI0);
        obj.ZAI=DAT.ZAI0(obj.burnIdx);
        obj.atomicMass=DAT.AMASS(obj.burnIdx);
        obj.decayEnergy=DAT.Q(obj.burnIdx);
        obj.halfLife=DAT.T12(obj.burnIdx);
        obj.oxState=valenceStates(obj.ZAI);
        [obj.ingToxicity,obj.inhToxicity]=ingAndInhTox(obj.ZAI);
        obj.hasNucData=hasXSData(obj.ZAI);
        obj.nuclideName=ZAI2Name(obj.ZAI);
        defDecMtx=sparse(DAT.decayMatrix);
        defDecMtx=defDecMtx(:,obj.burnIdx);
        defDecMtx=defDecMtx(obj.burnIdx,:);
        obj.defDecMtx=defDecMtx;
      else
        error('Nuclear data library not found.')
      end
    end
  end
  
end

