function [initialFrameArray,cineMetaData] = findBGTethered(exprPath,movNum,colorOffsets,windowParams)
%FINDBGTETHERED Summary of this function goes here
%   Detailed explanation goes here
LoadPhantomLibraries();
RegisterPhantom(true);

metaDataDir = dir(fullfile(exprPath,'*.cine'));
cineMetaData = getCinMetaData(fullfile(metaDataDir(1).folder,...
    metaDataDir(1).name));

% camNames = {'xy','yz','xz'};
camNames = {'yz','xz','xy'};
initialFrameArray = cell(1,3);
currImFrame = cell(1,3);
window_length = 60;
cineCenter = [round(cineMetaData.width/2),round(cineMetaData.height/2)];

movStrNum = sprintf('%03d',movNum);

bgCellFile = fullfile(fileparts(exprPath),'bg','bgCell.mat');
load(bgCellFile);

for camInd = 1:3
    % save background image
    % currBGDir = dir(fullfile(exprPath,'bg',[camNames{camInd},'*.cine']));
    % camFilenameBG = currBGDir(1).name;
    % % camFilenameBG = [camNames{camInd},'_001.cine'];
    % cinePath_bg = fullfile(exprPath,'bg',camFilenameBG);
    % cineDataBG = myOpenCinFile(cinePath_bg);
    % bgIm = myReadCinImage(cineDataBG,1);
    % myCloseCinFile(cineDataBG);
    bgIm = bgCell{camInd};

    % Load input im 
    currCamDir = dir(fullfile(exprPath,[camNames{camInd},'_',movStrNum,...
        '.cine']));
    cineFilename = currCamDir(1).name;
    cinePath = fullfile(exprPath,cineFilename);
    currCineData = myOpenCinFile(cinePath);
    currIm = myReadCinImage(currCineData,0);
    currImFrame{camInd} = currIm;
    
    tempIm = currIm;
    w_dims = windowParams(camInd,:);
    tempIm(w_dims(1):w_dims(2),w_dims(3):w_dims(4)) = ...
        bgIm(w_dims(1):w_dims(2),w_dims(3):w_dims(4))+colorOffsets(camInd);
    
    initialFrameArray{camInd} = tempIm;

    myCloseCinFile(currCineData);
end

UnregisterPhantom(); %Use this function when you finished your work
UnloadPhantomLibraries();
end



