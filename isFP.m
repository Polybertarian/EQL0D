function isfp = isFP(ZAIList)
%ISFP returns logicals saying is input ZAI corresponds to a fission product
lowThreshold=210000; highThreshold=730000;
FPList=ZAIList(((ZAIList>1000&ZAIList<3e4)|ZAIList>lowThreshold)&(ZAIList<highThreshold));
FPList(ismember(FPList,lowThreshold:1e4:highThreshold))=[];

isfp=ismember(ZAIList,FPList);

end
