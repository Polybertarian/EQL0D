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
    end
    
    methods
        function obj = NuclearObject()
            global DAT
            if exist('DAT','var')
                obj.ZAI=DAT.ZAI0;
                obj.atomicMass=DAT.AMASS;
                obj.decayEnergy=DAT.Q;
                obj.halfLife=DAT.T12;
                obj.oxState=valenceStates(obj.ZAI);
                [obj.ingToxicity,obj.inhToxicity]=ingAndInhTox(obj.ZAI);
                obj.hasNucData=hasXSData(obj.ZAI);
                obj.nuclideName=ZAI2Name(obj.ZAI);
            else
                error('Nuclear data library not found.')
            end
        end
    end
    
end

