% script that will run generateTetheredMovies over all cines for a
% particular genotype/type of experiment

fnExpMovNum = '(?<camName>[xyz]{2})_(?<movieNum>\d{3}).cine';
exprPath = 'Y:\Abby2\Feb7';
% allFolderNames = {'S01','S02','control haltereless','S01 haltereless','S02 haltereless'};
% allFolderPaths = fullfile(exprPath,allFolderNames,'fly*');
flyFolders = dir(fullfile(exprPath,'**','fly *'));

%% loop over folders and generate mp4s
for flyInd = 1:length(flyFolders)
    currFlyFolder = fullfile(flyFolders(flyInd).folder,flyFolders(flyInd).name);
    % load fly COM data
    try
        windowParams = load(fullfile(currFlyFolder,'bg_info','flyWindowParams.mat'));
        windowParams = windowParams.windowParams;
    catch
        continue
    end
    saveToPath = fullfile(currFlyFolder,'mp4');

    % get movie numbers
    cineDir = dir(fullfile(currFlyFolder,'*.cine'));
    renamed_out = regexp({cineDir.name},fnExpMovNum,'names');
    renamedInds = ~(cellfun(@isempty,renamed_out));
    movsRenamed = renamed_out(renamedInds);
    renamedStruct = [movsRenamed{:}];
    movNumsList = str2double({renamedStruct.movieNum});

    % loop over each movie and generate mp4
    for movNum = unique(movNumsList)
        % check for triplet
        tripletCheck = sum(movNumsList==movNum);
        if tripletCheck ~= 3
            continue
        end
        tic
        generateTetheredMovies(currFlyFolder,windowParams,saveToPath,movNum)
        toc
    end
end
