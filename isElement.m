function [isElementList,oldIdx] = isElement(elementList,zaiList)
%ISELEMENT returns logical if some isotopes of ZAILIST are elements of
%ELEMENTSLIST

isElementList=false(size(zaiList));
oldIdx=zeros(size(zaiList));
elementList=unique(elementList);

for i=1:numel(elementList)
    idx=zaiList>elementList(i)*1e4&zaiList<(elementList(i)+1)*1e4;
    isElementList(idx)=true;
    oldIdx(idx)=i;
end
oldIdx(oldIdx==0)=[];

return
end
