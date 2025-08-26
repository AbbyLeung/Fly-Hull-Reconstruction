%{
Function that converts 4 cine files to an mp4. Will throw an error if the
cine params are different.

Inputs
cineFilenames: cell array where each element is the full path of each cine 
file. Should be in order, so the first path is cam1, second path is cam2,
and so on.
saveToPath: path the mp4 will be saved to. Ideally, this will be the mp4
folder.
movNum: movie number that identifies which movie the mp4 corresponds to.
This is used for naming the mp4 file. 
%}

function [] = generateTetheredMovies(exprPath,windowParams,saveToPath,movNum)
% set up paths
camNames = {'yz','xz','xy'};

cineFilenames = strcat(camNames,'_',sprintf('%03d',movNum),'.cine');
cinePaths = fullfile(exprPath,cineFilenames);
mp4Filename = ['Movie_',int2str(movNum),'.mp4'];
mp4Path = fullfile(saveToPath,mp4Filename);
numCams = length(cinePaths);

LoadPhantomLibraries();
RegisterPhantom(true); 

% if times and frame size are different, throw error
metaDataCell = cell(1,length(cinePaths));
keepFields = {'firstImage','lastImage','width','height','framerate'};
for cineInd = 1:length(cinePaths)
    currMetaData = getCinMetaData(cinePaths{cineInd});
    currFieldNames = fieldnames(currMetaData);
    removeFields = currFieldNames(~ismember(currFieldNames,keepFields));
    metaDataCell{cineInd} = rmfield(currMetaData,removeFields);
end

if ~isequal(metaDataCell{:})
    error('cine parameters are not the same across all videos.')
end
clear metaDataCell
% get metadata and set video params
metaData = getCinMetaData(cinePaths{1});
framerate = double(metaData.frameRate);
firstImNo = metaData.firstImage;
lastImNo = metaData.lastImage;
frameNoArray = firstImNo:lastImNo;
timeArray = floor((frameNoArray/framerate)*1000);
N_frames = length(timeArray);
startFrameNum = firstImNo;

imWidth = 200;
imHeight = 200;

% set text for time stamps
timeVals = (firstImNo:lastImNo) * (1/framerate) * 1000;
timeVals = floor(timeVals);
timeValsText = cell(1,length(timeVals));

for timeInd = 1:length(timeVals)
    timeValsText{timeInd}=[int2str(timeVals(timeInd)),'ms'];
end

% Get cine data for opening the cine files
cineDataCell = cell(1,numCams);
for cineDataInd = 1:numCams
    cineDataCell{cineDataInd} = myOpenCinFile(cinePaths{cineDataInd});
end

writerObj=VideoWriter(mp4Path,'MPEG-4');
writerObj.Quality = 100;
open(writerObj);
allWidthInds = {1:imWidth,imWidth+1:imWidth*2,imWidth*2+1:imWidth*3};
allHeightInds = {1:imHeight,1:imHeight,1:imHeight};

% write video in chunks since vid variable is large
chunkSize = 365;
N_chunks = ceil(N_frames/chunkSize);
frameInd = 1;
for currChunk = 1:N_chunks
    clear vid
    if (currChunk==N_chunks)
        currChunkSize = N_frames - chunkSize*(currChunk-1);
    else
        currChunkSize = chunkSize;
    end
    % define vid, which is the image stack for the movie in this chunk
    vid.height = imHeight;
    vid.width = imWidth*3;
    fr = zeros(vid.height, vid.width, 3, 'uint8');
    vid.frames(1:currChunk) = struct('cdata',fr);

        for camInd = 1:numCams
            vidInd = 0;
            % loop through frames in chunk
            for currFrameNum = startFrameNum:(startFrameNum+currChunkSize-1)
                vidInd = vidInd+1;
                currIm = myReadCinImage(cineDataCell{camInd},currFrameNum);
                COM = floor(mean([windowParams(camInd,1:2);windowParams(camInd,3:4)],2));
                heightInds = allHeightInds{camInd};
                widthInds = allWidthInds{camInd};
                pad = [imWidth/2,imHeight/2];

                currFrame = currIm((COM(1)-pad(1)):(COM(1)+pad(1)-1),...
                    (COM(2)-pad(2)):(COM(2)+pad(2)-1));

                % imbed timer in cam1 view
                if camInd == 1
                    currFrame = insertText(currFrame,[1,1],...
                        timeValsText{frameInd},BoxOpacity=1,FontSize=14);
                    currFrame = currFrame(:,:,1);
                    frameInd = frameInd+1;
                end
                vid.frames(vidInd).cdata(heightInds,widthInds,1) = im2uint8(currFrame);
                vid.frames(vidInd).cdata(heightInds,widthInds,2) = im2uint8(currFrame);
                vid.frames(vidInd).cdata(heightInds,widthInds,3) = im2uint8(currFrame);
            end
        end
    for writeFrameInd = 1:vidInd
        writeVideo(writerObj,vid.frames(writeFrameInd).cdata);
    end
    startFrameNum = currFrameNum+1;
end

close(writerObj);
UnregisterPhantom();
UnloadPhantomLibraries();
end

