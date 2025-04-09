% script to load data into one data structure
exprFolder = 'Y:\Abby2\Feb7\Control haltereless';
flyLine = 'Control haltereless';
flyDir = dir(fullfile(exprFolder,'fly*'));
folderNames = {flyDir.name};
datastruct = struct();
structInd = 1;

% loop through each folder and trial to extract data
for flyInd = 1:length(folderNames)
    currFlyStr = regexp(folderNames{flyInd},'\d*','Match');
    currFlyNum = str2double(currFlyStr);

    analysisFolder = fullfile(exprFolder,folderNames{flyInd},'Analysis');
    trialsDir = dir(fullfile(analysisFolder,'Expr*'));

    if isempty(trialsDir)
        continue
    end
    
    for trialInd = 1:length(trialsDir)
       trialFolderName = trialsDir(trialInd).name;
       currTrialStr = regexp(trialFolderName,'\d*','Match'); % get num from filename
       currTrialNum = str2double(currTrialStr{end});
       % add info to data struct
       datastruct(structInd).flyLine = flyLine;
       datastruct(structInd).flyNumber = currFlyNum;
       datastruct(structInd).trialNum = currTrialNum;
        
       % load analysis file
       currFolder = fullfile(trialsDir(trialInd).folder,trialFolderName);
       mcDir = dir(fullfile(currFolder,'*manually_corrected.mat'));

       try 
           if ~isempty(mcDir)
               dataFilename = mcDir.name;
               currData = load(fullfile(currFolder,dataFilename));
               currData = currData.data;
           else % assume the data is cleaned already (part of reconstruction)
                dataFilename = [trialFolderName,'_cleaned.mat'];
                currData = load(fullfile(currFolder,dataFilename));
                currData = currData.data_cleaned;
           end
       catch errorMsg
            disp(errorMsg)
            continue
       end

       % currData = currData.data_cleaned;
       bodyAngles = currData.anglesLabFrameSmooth(:,[1,2,9]); %yaw, pitch, roll
       wingAngles = currData.anglesBodyFrameSmooth;
       wingAngles(:,3) = wingAngles(:,3)*(-1);

       wingAnglesR = wingAngles(:,[3,4,5]);
       wingAnglesL = wingAngles(:,[6,7,8]);

       % get current time in seconds
       firstFrame = currData.params.startTrackingTime;
       lastFrame = currData.params.endTrackingTime;
       fps = currData.params.fps;

       time = (firstFrame:lastFrame)*(1/fps);


       % get fly shadow
       % tic
       % flyShadow = getXYview(currData);
       % toc

       datastruct(structInd).t  = time;
       datastruct(structInd).bodyYaw = bodyAngles(:,1);
       datastruct(structInd).bodyPitch = bodyAngles(:,2);
       datastruct(structInd).bodyRoll = bodyAngles(:,3);
       datastruct(structInd).wing_smoothR = wingAnglesR';
       datastruct(structInd).wing_smoothL = wingAnglesL';
       % datastruct(structInd).flyShadow = flyShadow;

       clear currData
       structInd = structInd + 1;
    end  
end

S01Datastruct = datastruct;
save('S01_datastruct.mat','S01Datastruct')
