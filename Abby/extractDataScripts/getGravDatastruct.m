% script to load data into one data structure
exprFolder = 'D:\Gravity Sensing\Fly 01\Intact_Light';
flyLine = 'Intact In Light Wildtype';
flyDir = dir(fullfile(exprFolder,'fly*'));
% folderNames = {flyDir.name};
datastruct = struct();
structInd = 1;

% loop through each folder and trial
for flyInd = 1%:length(folderNames)
    % currFlyStr = regexp(folderNames{flyInd},'\d*','Match');
    % currFlyNum = str2double(currFlyStr);
    currFlyNum = 1;

    % analysisFolder = fullfile(exprFolder,folderNames{flyInd},'Analysis');
    analysisFolder = fullfile(exprFolder,'Analysis');
    trialsDir = dir(fullfile(analysisFolder,'Expr*'));

    if isempty(trialsDir)
        continue
    end
    
    for trialInd = 1:length(trialsDir)
       trialFolderName = trialsDir(trialInd).name;
       currTrialStr = regexp(trialFolderName,'\d*','Match'); % get num from filename
       currTrialNum = str2double(currTrialStr{end});
       % load analysis file
       currFolder = fullfile(trialsDir(trialInd).folder,trialFolderName);

       mcDir      = dir(fullfile(currFolder,'*manually_corrected.mat')) ;
       cleanedDir = dir(fullfile(currFolder,'*cleaned.mat')) ;
       testDir    = dir(fullfile(currFolder,'*test.mat')) ;

       if ~isempty(mcDir)
           dataFilename = mcDir.name ;
           varName      = 'data' ;
           manualCorr   = true ;
       elseif ~isempty(cleanedDir)
           dataFilename = cleanedDir.name ;
           varName      = 'data_cleaned' ;
           manualCorr   = false ;
       elseif ~isempty(testDir)
           dataFilename = testDir.name ;
           varName      = 'data' ;
           manualCorr   = false ;
       else
           fprintf('No data file found in %s, skipping\n', trialFolderName) ;
           continue
       end

       try
           raw      = load(fullfile(currFolder, dataFilename)) ;
           currData = raw.(varName) ;
       catch errorMsg
           disp(errorMsg)
           continue
       end

       % add info to data struct (only after confirming data loaded)
       datastruct(structInd).flyLine         = flyLine ;
       datastruct(structInd).flyNumber       = currFlyNum ;
       datastruct(structInd).trialNum        = currTrialNum ;
       datastruct(structInd).dataPath        = currFolder ;
       datastruct(structInd).ManualCorr      = manualCorr ;
       datastruct(structInd).ManualCorrRange = [0, 30] ;

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
% 
ctrlGrav1_light = datastruct;
save('ctrlGrav1_light.mat','ctrlGrav1_light')
