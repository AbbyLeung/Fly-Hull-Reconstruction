function [initialFrameArray,cineMetaData] = findBGTethered(exprPath)
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
% offset_vec = [13,13,15];
offset_vec = [0,0,0];
window_length = 100;
w_dims = [round(cineMetaData.width/2)-window_length,round(cineMetaData.height/2)+window_length];

% this offset vec is because I took the background image on a different day
% from the experiment because I didn't know better. I shouldn't have to do
% this in the future though.

for camInd = 1:3
    % save background image
    currBGDir = dir(fullfile(exprPath,'bg',[camNames{camInd},'*.cine']));
    camFilenameBG = currBGDir(1).name;
    % camFilenameBG = [camNames{camInd},'_001.cine'];
    cinePath_bg = fullfile(exprPath,'bg',camFilenameBG);
    cineDataBG = myOpenCinFile(cinePath_bg);
    bgIm = myReadCinImage(cineDataBG,1);
    myCloseCinFile(cineDataBG);

    % Load input im 
    currCamDir = dir(fullfile(exprPath,[camNames{camInd},'*.cine*']));
    cineFilename = currCamDir(1).name;
    cinePath = fullfile(exprPath,cineFilename);
    currCineData = myOpenCinFile(cinePath);
    currIm = myReadCinImage(currCineData,0);
    currImFrame{camInd} = currIm;
    
    tempIm = currIm;

    tempIm(w_dims(1):w_dims(2),w_dims(1):w_dims(2)) = bgIm(w_dims(1):w_dims(2),w_dims(1):w_dims(2))-offset_vec(camInd);
    initialFrameArray{camInd} = tempIm;

    myCloseCinFile(currCineData);
end

UnregisterPhantom(); %Use this function when you finished your work
UnloadPhantomLibraries();
end

