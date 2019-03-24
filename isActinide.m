function isAct = isActinide (NuclideZAI)
%ISACTINIDE tests if target nuclide is an actinide or not
ZAIThreshold=88e4;
isAct=NuclideZAI>ZAIThreshold|(NuclideZAI<112&NuclideZAI>ZAIThreshold/1e4);
end
