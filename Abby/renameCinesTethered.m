% Script to rename cine files from the trigger time to a number. This
% also checks if there are cines already renamed. 
falseTriggerFlag = false;
errorFlag = false;

fnExp1 = ['(?<camName>[xyz]{2})_Y(?<year>\d{4})',...
    '(?<month>\d{2})(?<day>\d{2})H(?<hour>\d{2})(?<minute>\d{2})', ...
    '(?<second>\d+.\d+)']; % two digits in day

fnExp3 = ['(?<camName>[xyz]{2})_Y(?<year>\d{4})',...
    '(?<month>\d{2}) (?<day>\d{1})H(?<hour>\d{2})(?<minute>\d{2})', ...
    '(?<second>\d+.\d+)']; % one digit in day AND ALSO a space between the month and day
fnExp = [fnExp1,'|',fnExp3]; 

% reg expression for movies that are already renamed
fnExpMovNum = '(?<camName>[xyz]{2})_(?<movieNum>\d{3}).cine';
pathToWatch =  'Y:\Abby\tethered_3cam_data\02_20012025\';

% get expr num
currFilename = split(pathToWatch,'\');
exprNum = split(currFilename(end),'_');
exprNum = str2double(exprNum{1});
exprNumStr = num2str(exprNum,'%03.f');

% make folders
mp4Path = fullfile(pathToWatch,'mp4');
if ~exist(mp4Path,'dir')
    mkdir(mp4Path)
end

unsortedPath = fullfile(pathToWatch,'Analysis','Unsorted');
if ~exist(unsortedPath,'dir')
    mkdir(unsortedPath)
end

falseTriggerPath = fullfile(pathToWatch,'Possible False Triggers');
if ~exist(falseTriggerPath,'dir')
    mkdir(falseTriggerPath);
end

cineDirs = cell(1,3);
datetimes = cell(1,3);

camNames = {'xy','yz','xz'};

% get trigger times from each camera
for camInd = 1:3
    camFilename = [camNames{camInd},'*.cine'];
    currDir = dir(fullfile(pathToWatch,camFilename));
    filenameLengths = cellfun(@length,{currDir.name});
    renameInds = filenameLengths > 20; % short filenames are already renamed
    cineDirs{camInd} = currDir(renameInds);

    filenames = {cineDirs{camInd}.name};

    for fileInd = 1:length(cineDirs{camInd})
        out_names = regexp(filenames,fnExp,'names') ;
        currName = out_names{fileInd};

        if length(currName.day) == 1
            currName.day = ['0',currName.day];
        end

        datetime_str = [currName.year,currName.month,currName.day,...
            currName.hour,currName.minute,currName.second(1:7)];
        cineDirs{camInd}(fileInd).triggerTime = datetime(datetime_str,'InputFormat',...
            'yyyyMMddHHmmss.SSSS');
        cineDirs{camInd}(fileInd).camName = currName.camName;

        datetimes{camInd} = [cineDirs{camInd}.triggerTime];
    end
end

% get current movie number
allMovDir = dir(fullfile(pathToWatch,'*.cine'));
movNum_outputs = regexp({allMovDir.name},fnExpMovNum,'names');
movNumInds = ~(cellfun(@isempty,movNum_outputs));

if sum(movNumInds) == 0
    movNum = 0;
else
    movNumsRenamed = movNum_outputs(movNumInds);
    movNumStruct = [movNumsRenamed{:}]; % restructure into superior struct array
    movNumsList = str2double({movNumStruct.movieNum});
    
    movNum = max(movNumsList);
end

% compare trigger times across the 3 cams
for ind1 = 1:length(datetimes{1})
    [timeDiffs_2,ind2] = min(abs((datetimes{1}(ind1)-datetimes{2})));
    [timeDiffs_3,ind3] = min(abs((datetimes{1}(ind1)-datetimes{3})));

    tol = duration(0,0,0.5); % 0.5 second tolerance
    if timeDiffs_2 < tol && timeDiffs_3 < tol
        movNum = movNum + 1;
        movNumStr = num2str(movNum,'%03.f');

        origFilenames = {cineDirs{1}(ind1).name,cineDirs{2}(ind2).name,...
            cineDirs{3}(ind3).name};
        origPaths = fullfile(pathToWatch,origFilenames);
        renameFilenames = strcat(camNames,['_',movNumStr,'.cine']);
        cinePaths = fullfile(pathToWatch,renameFilenames);

        % rename files
        movefile(origPaths{1},cinePaths{1}); movefile(origPaths{2},cinePaths{2});
        movefile(origPaths{3},cinePaths{3});

        % also rename the .xml files
        [~,cam1Name,~] = fileparts(origPaths{1});
        [~,cam2Name,~] = fileparts(origPaths{2});
        [~,cam3Name,~] = fileparts(origPaths{3});
        
        xmlFilenames = strcat(camNames,['_',movNumStr,'.xml']);
        xmlPaths = fullfile(pathToWatch,xmlFilenames);

        movefile(fullfile(pathToWatch,[cam1Name,'.xml']),xmlPaths{1});
        movefile(fullfile(pathToWatch,[cam2Name,'.xml']),xmlPaths{2});
        movefile(fullfile(pathToWatch,[cam3Name,'.xml']),xmlPaths{3});
    end
end

