function isfp = isFP(ZAIList)
%ISFP returns logicals saying is input ZAI corresponds to a fission product
lowThreshold=21e4;
highThreshold=73e4;
FPList=ZAIList(ZAIList<3e4|ZAIList>lowThreshold&ZAIList<highThreshold);
FPList(ismember(FPList,lowThreshold:1e4:highThreshold))=[];

isfp=ismember(ZAIList,FPList);

end

