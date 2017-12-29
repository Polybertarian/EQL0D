classdef Rep
    %REP Reprocessing stream in EQL0D
    
    properties
        name    
        srcMat   
        dstMat  
        elements
        share
        frequency
        rate
        type    
        mode    
    end
    
    methods
        function obj = Rep(name,srcMat,dstMat,elements,share,rate,type)
            global DAT
            obj.name=name;
            obj.srcMat=srcMat;
            obj.dstMat=dstMat;
            if(isempty(share))
                obj.share=ones(size(elements));
            elseif(length(share)~=length(elements))
                error('Error: share and elements vector length mismatch.');
            else
                obj.share=share;
            end
            if(~iscolumn(obj.share))
                obj.share=obj.share';
            end
            if(~iscolumn(obj.elements))
                obj.elements=obj.elements';
            end
            if(any(obj.elements<111))
                [idx1,idx2]=isElement(obj.elements,DAT.ZAI0);
                obj.elements=DAT.ZAI0(idx1);
                obj.share=obj.share(idx2);
            end
            obj.elements=elements;
            obj.type=type;
            obj.frequency=1;
            if(ischar(rate))
                obj.rate=1;
                if(ismember(rate,{'keep','keepAM','keepAFPM','keepTotM','keepAA','keepAFPM','keepTotA'}))
                    obj.mode=rate;
                else
                    error(['Reprocessing stream mode ' rate ' of stream ' obj.name ' not recognized!']);
                end
            elseif(isnumeric(rate))
                obj.rate=rate;
                obj.mode='remove';
            end
        end
        function srcIdx = findSrc(obj,matList)
            if(isempty(obj.srcMat))
                srcIdx=0;
            else
                srcIdx=find(ismember(matList,obj.srcMat));
            end
        end
        function dstIdx = findDst(obj,matList)
            if(isempty(obj.dstMat))
                dstIdx=0;
            else
                dstIdx=find(ismember(matList,obj.dstMat));
            end
        end
        function obj = adaptElements(obj,ZAIList)
            if(any(obj.elements<111))
                [idx1,idx2]=isElement(obj.elements,ZAIList);
                obj.elements=ZAIList(idx1);
                obj.share=obj.share(idx2);
            end
        end
    end
    
end

