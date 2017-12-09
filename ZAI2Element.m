function elementNumber = ZAI2Element(nuclideZAI)
%ZAI2ELEMENT gives the corresponding element number for a given ZAI
strArray=cellstr(strtrim(num2str(nuclideZAI)));
elementNumber=zeros(size(nuclideZAI));
for i=1:length(strArray)
    elementNumber(i)=str2double(strArray{i}(1:end-4));
end
end

