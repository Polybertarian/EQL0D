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
  end
  
  methods
    function obj = NuclearObject()
      global DAT
      if exist('DAT','var')
        idx=isProduced(DAT.libraryName,DAT.ZAI0);
        obj.ZAI=DAT.ZAI0(idx);
        obj.atomicMass=DAT.AMASS(idx);
        obj.decayEnergy=DAT.Q(idx);
        obj.halfLife=DAT.T12(idx);
        obj.oxState=valenceStates(obj.ZAI);
        [obj.ingToxicity,obj.inhToxicity]=ingAndInhTox(obj.ZAI);
        obj.hasNucData=hasXSData(obj.ZAI);
        obj.nuclideName=ZAI2Name(obj.ZAI);
        defDecMtx=sparse(DAT.decayMatrix);
        defDecMtx=defDecMtx(:,idx);
        defDecMtx=defDecMtx(idx,:);
        obj.defDecMtx=defDecMtx;
      else
        error('Nuclear data library not found.')
      end
    end
  end
  
end

