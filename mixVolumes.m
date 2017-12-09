function newDens = mixVolumes(Mat,Rep,Volume)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if(Volume==0)
    newDens=Mat.atDens;
elseif(Volume<0&abs(Volume)>Mat.volume)
    newDens=1e-50*ones(length(Mat.ZAI),1);
else
    N=zeros(length(Mat.ZAI),1);
    N(ismember(Mat.ZAI,Rep.elements))=Rep.share;
    newDens=Mat.atDens/(1+Volume/Mat.volume)+N/(1+Mat.volume/Volume);
    newDens(newDens<0)=0.0;
end

return
end