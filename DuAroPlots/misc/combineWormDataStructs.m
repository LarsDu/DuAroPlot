%Combine wormData structs
function [combinedWormData] = combineWormDataStructs( Data1, Data2)
combinedWormData = struct;
combinedWormData.spotNum = [Data1.myWormDataStruct.spotNum;Data2.myWormDataStruct.spotNum];
combinedWormData.U = [Data1.myWormDataStruct.U;Data2.myWormDataStruct.U];
combinedWormData.L = [Data1.myWormDataStruct.L;Data2.myWormDataStruct.L];
