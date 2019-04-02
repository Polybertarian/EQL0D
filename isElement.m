function [isElementList,oldIdx] = isElement(elementList,zaiList)
    %ISELEMENT returns logical if some isotopes of ZAILIST are elements of
    %ELEMENTSLIST

    isElementList=false(size(zaiList)); oldIdx=[];

    for i=1:numel(elementList)
        if ismember(elementList(i),zaiList)
            idx=zaiList==elementList(i);
        else
            idx=zaiList>elementList(i)*1e4&zaiList<(elementList(i)+1)*1e4;
        end
        isElementList(idx)=true;
        oldIdx(end+1:end+numel(find(idx)))=i;
    end
end
