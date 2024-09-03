function [initialFrameArray] = findBGTethered(exprPath)
%FINDBGTETHERED Summary of this function goes here
%   Detailed explanation goes here
LoadPhantomLibraries();
RegisterPhantom(true);

% exprPath = 'C:\Users\Abby\Cornell\Abby_test_tethered_code\01_23082024';

% cinePath = fullfile(exprPath,'xy_001.cine');
% cineData = myOpenCinFile(cinePath);
% cineMetaData = getCinMetaData(cinePath);

% camNames = {'xy','yz','xz'};
camNames = {'yz','xz','xy'};
initialFrameArray = cell(1,3);
currImFrame = cell(1,3);
offset_vec = [13,13,15];

for camInd = 1:3
    % save background image
    camFilenameBG = [camNames{camInd},'_3.cine'];
    cinePath_bg = fullfile(exprPath,'bg',camFilenameBG);
    cineDataBG = myOpenCinFile(cinePath_bg);
    bgIm = myReadCinImage(cineDataBG,1);
    myCloseCinFile(cineDataBG);

    % Load input im 
    cineFilename = [camNames{camInd},'_001.cine'];
    cinePath = fullfile(exprPath,cineFilename);
    currCineData = myOpenCinFile(cinePath);
    currIm = myReadCinImage(currCineData,0);
    currImFrame{camInd} = currIm;
    
    tempIm = currIm;
    tempIm(300:450,250:400) = bgIm(300:450,250:400)-offset_vec(camInd);
    initialFrameArray{camInd} = tempIm;

    myCloseCinFile(currCineData);
end

UnregisterPhantom(); %Use this function when you finished your work
UnloadPhantomLibraries();
end

