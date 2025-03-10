function [flyShadowSum] = getXYview(data)
%GETXYVIEW Summary of this function goes here
%   Detailed explanation goes here

res = data.res;
frameNumList = res(:,1);
frameNums = unique(frameNumList);
flyShadowSum = zeros(1,length(frameNums));

for frameInd = 1:length(frameNums)
    currpts = res(frameNumList==frameNums(frameInd),2:3);
    numPts = unique(currpts,'rows');
    flyShadowSum(frameInd) = height(numPts);
end

