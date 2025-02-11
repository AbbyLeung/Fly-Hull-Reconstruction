%% first get background images (only need to do this one time)
bgFolder = 'Y:\Abby2\Feb7\bg';
camNames = {'yz','xz','xy'};
camFilenames = strcat(camNames,'_1.cine');
camPaths = fullfile(bgFolder,camFilenames);

% check if bg cell is already created
bgCellDir = dir(fullfile(bgFolder,'bgCell*'));

if ~isempty(bgCellDir) % bg cell already exists
    load(fullfile(bgFolder,bgCellDir(1).name))
else
    LoadPhantomLibraries();
    RegisterPhantom(true);

    bgCell = cell(1,3);

    for camInd = 1:3
        currCineData = myOpenCinFile(camPaths{camInd});
        currIm = myReadCinImage(currCineData,0);
        bgCell{camInd} = currIm;
        myCloseCinFile(currCineData);
    end

    UnregisterPhantom(); %Use this function when you finished your work
    UnloadPhantomLibraries();

    save(fullfile(bgFolder,'bgCell.mat'),'bgCell')
end

% make mask
dims = size(bgCell{1});
mask = false(dims);
mask(100:200,650:700) = true;
colorOffsets = zeros(1,3);

% get one frame from each camera view to adjust color in background for
% each fly
% bgInfoPath = 'Y:\Abby2\Feb7\S01\fly 07\bg_info';
load(fullfile(bgInfoPath,'bgParams.mat'))
for camInd = 1:3
    avgBG = mean(bgCell{camInd}(mask));
    avgIm = mean(imgCell{camInd}(mask));
    colorOffsets(camInd) = floor(avgIm-avgBG);
end

%% get frames across one wingbeat
LoadPhantomLibraries();
RegisterPhantom(true);

% cineFilenames = strcat(camNames,'_001.cine');
% cinePath = fullfile('Y:\Abby2\Feb7\S01\fly 10',cineFilenames);
cinePath = cineFilenames;
frames = struct();
windowParams = zeros(3,4);

for camInd = 1:3
    currCineData = myOpenCinFile(cinePath{camInd});
    
    for ind = 1:36
        currIm = myReadCinImage(currCineData,ind-1);
        frames(ind).(camNames{camInd}) = currIm;
    end
    
    myCloseCinFile(currCineData);
    
    framesSum = imbinarize(bgCell{camInd}+colorOffsets(camInd)-frames(1).(camNames{camInd}));
    for i = 2:36
        framesSum = framesSum+(imbinarize(bgCell{camInd}+colorOffsets(camInd)-frames(i).(camNames{camInd})));
    end
    
    framesSum(framesSum==36)=0;
    BW = imbinarize(framesSum);
    BW2 = imfill(BW,'holes');
    % imshow(BW2)
    
    % find the largest 2 connected components, which will be the wings
    CC = bwconncomp(BW2);
    compSize = cellfun(@length, CC.PixelIdxList);
    [test,compInds] = sort(compSize);
    wingInds = [compInds(end-1),compInds(end)];

    wingLinearInds = [CC.PixelIdxList{wingInds(1)};CC.PixelIdxList{wingInds(2)}];
    [row,col] = ind2sub(dims,wingLinearInds);
    rowMax = max(row) + 5;
    rowMin = min(row) - 5;
    colMax = max(col) + 5;
    colMin = min(col) - 5;

    windowParams(camInd,:) = [rowMin,rowMax,colMin,colMax];
end
UnregisterPhantom(); %Use this function when you finished your work
UnloadPhantomLibraries();

save(fullfile(bgInfoPath,'flyWindowParams.mat'),'colorOffsets','windowParams')

% testFrame = frames(1).XZ;
% testFrame(rowMin:rowMax,colMin:colMax) = 256;
% 
% 
% currCamInd = 2;
% for i = 1:36
%     pause(0.5)
%     imshow(frames(i).(camNames{currCamInd})(windowParams(currCamInd,1):windowParams(currCamInd,2),...
%         windowParams(currCamInd,3):windowParams(currCamInd,4)))
%     impixelinfo
% end
% 
% 








