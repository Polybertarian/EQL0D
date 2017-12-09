classdef Rep
    %REP Reprocessing stream in EQL0D
    
    properties
        name
        srcMat
        dstMat
        elements
        share
        rate
        type
        mode
    end
    
    methods
        function obj = Rep(name,srcMat,dstMat,elements,share,rate,type,mode)
           obj.name=name;
           obj.srcMat=srcMat;
           obj.dstMat=dstMat;
           obj.elements=elements;
           obj.share=share;
           obj.rate=rate;
           obj.type=type;
           obj.mode=mode;
        end
    end
    
end

