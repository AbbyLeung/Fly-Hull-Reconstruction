function [initialFrameArray,cineMetaData] = findBGTethered(exprPath,movNum)
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
window_length = 100;
cineCenter = [round(cineMetaData.width/2),round(cineMetaData.height/2)];
dims = cell(1,3);

dims{1} = cineCenter + [-window_length,window_length];
dims{2} = cineCenter + [-window_length,window_length];
dims{3} = cineCenter + [-window_length,window_length];

% this offset vec is because I took the background image on a different day
% from the experiment because I didn't know better. I shouldn't have to do
% this in the future though.

movStrNum = sprintf('%03d',movNum);

colorCorr = [-3,5,5];
widthOffset = [0,-150,0];
heightOffset = [0,0,-70];
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
    currCamDir = dir(fullfile(exprPath,[camNames{camInd},'_',movStrNum,...
        '.cine']));
    cineFilename = currCamDir(1).name;
    cinePath = fullfile(exprPath,cineFilename);
    currCineData = myOpenCinFile(cinePath);
    currIm = myReadCinImage(currCineData,0);
    currImFrame{camInd} = currIm;
    
    w_dims = dims{camInd};

%     if camInd == 1 || camInd == 2
%         tempIm = currIm;
%         widthOffset = -150; heightOffset = 0;
%         tempIm((w_dims(1):w_dims(2))+heightOffset,...
%             (w_dims(1):w_dims(2))+widthOffset) = bgIm((w_dims(1):w_dims(2))+heightOffset,...
%             (w_dims(1):w_dims(2))+widthOffset)+5;
%         imshow(tempIm)
%     else
%         tempIm = currIm;
%         tempIm(w_dims(1):w_dims(2),w_dims(1):w_dims(2)) = bgIm(w_dims(1):w_dims(2),w_dims(1):w_dims(2))-3;
%         imshow(tempIm);
%     end
%     

    tempIm = currIm;
    tempIm((w_dims(1):w_dims(2))+heightOffset(camInd),...
        (w_dims(1):w_dims(2))+widthOffset(camInd)) = bgIm((w_dims(1):w_dims(2))+heightOffset(camInd),...
        (w_dims(1):w_dims(2))+widthOffset(camInd))+colorCorr(camInd);
    
    initialFrameArray{camInd} = tempIm;

    myCloseCinFile(currCineData);
end

UnregisterPhantom(); %Use this function when you finished your work
UnloadPhantomLibraries();
end

