% Script to run analysis on renamed files. Need to modify the cnie2sparse
% code to make everything work.

pathToWatch = 'C:\Users\Abby\Cornell\Abby_test_tethered_code\01_23082024\' ; 
pathStruct = generatePathStruct(pathToWatch) ;
ExprNum = pathStruct.ExprNum;

clustFlag = true ; % shich version of analysis script to run
largePertFlag = false  ; % is it a large perturbation?
removeLegsFlag = true ; % try to remove legs in binary threshold?
alignBBoxFlag = false ; % try to align images to avoid clipping?

movNum = 1;
% 
% flyAnalysisMain(movNum, ExprNum, pathStruct, clustFlag, ...
%                 largePertFlag, removeLegsFlag) ;

bgCellArray = cell(1,3);

cinePath = fullfile('C:\Users\Abby\Cornell\Abby_test_tethered_code\01_23082024','xy_001.cine');
cineData = myOpenCinFile(cinePath);
cineMetaData = getCinMetaData(cinePath);

inputIm=myReadCinImage(cineData, -500);
myCloseCinFile(cineData);
clear cineData
clear cineMetaData

camNames = {'xy','yz','xz'};

for camInd = 1:3
    camFilename = [camNames{camInd},'_3.cine'];
    cinePath_bg = fullfile('C:\Users\Abby\Cornell\Abby_test_tethered_code\01_23082024','bg',camFilename);
    cineData = myOpenCinFile(cinePath_bg);
    bgIm = myReadCinImage(cineData,1);

    bgCellArray{camInd} = bgIm;
    myCloseCinFile(cineData);
end


%% xy params
% bg_color = 125;
% im1 = inputIm;
% im1(1:350,:) = bg_color;
% im1(450:end,:) = bg_color;
% im1(:,1:250) = bg_color;
% im1(:,400:end) = bg_color;
% 
% % imshow(im1<100)
% 
% mask=bwareaopen(im1<100,50);

%% xz params
% bg_color = 125;
% im1 = inputIm;
% im1(1:300,:) = bg_color;
% im1(450:end,:) = bg_color;
% im1(:,1:250) = bg_color;
% im1(:,400:end) = bg_color;
% 
% % imshow(im1)
% 
% % imshow(im1<110)
% 
% mask=bwareaopen(im1<100,50);
% imshow(mask)

%% yz params
% bg_color = 125;
% im1 = inputIm;
% im1(1:300,:) = bg_color;
% im1(450:end,:) = bg_color;
% im1(:,1:250) = bg_color;
% im1(:,400:end) = bg_color;

% imshow(im1)

% imshow(im1<110)
% 
% mask=bwareaopen(im1<100,50);
% imshow(mask)



