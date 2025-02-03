% Script to run analysis on renamed files. Need to modify the cnie2sparse
% code to make everything work.

%% Setup up paths for analysis folder
pathToWatch = 'Y:\Abby2\02_20012025\' ; 
pathStruct = generatePathStruct(pathToWatch) ;
ExprNum = pathStruct.ExprNum;
camNamesList = {'xy','xz','yz'};

%% rename matching cines
cineDir = dir(fullfile(pathToWatch,'*.cine'));

fnExp1 = ['(?<camName>[xyz]{2})_Y(?<year>\d{4})',...
    '(?<month>\d{2})(?<day>\d{2})H(?<hour>\d{2})(?<minute>\d{2})', ...
    '(?<second>\d+.\d+)']; % two digits in day

fnExp2 = ['(?<camName>[xyz]{2})_Y(?<year>\d{4})',...
    '(?<month>\d{2}) (?<day>\d{1})H(?<hour>\d{2})(?<minute>\d{2})', ...
    '(?<second>\d+.\d+)']; % one digit in day AND ALSO a space between the month and day

fnExpMovNum = '(?<camName>[xyz]{2})_(?<movieNum>\d{3}).cine';
fnExp = [fnExp1,'|',fnExp2]; 

out_names = regexp({cineDir.name},fnExp,'names');
renamingInd = ~cellfun(@isempty,out_names);
original_filenames = {cineDir(renamingInd).name};
% restructure into struct array
names_struct = [out_names{:}];

for i = 1:length(names_struct)
    names_struct(i).filename = original_filenames{i};
end

% get current movie number
renamed_out = regexp({cineDir.name},fnExpMovNum,'names');
renamedInds = ~(cellfun(@isempty,renamed_out));

if sum(renamedInds) == 0
    movNum = 0;
else
    movsRenamed = renamed_out(renamedInds);
    renamedStruct = [movsRenamed{:}];
    movNumsList = str2double({renamedStruct.movieNum});

    movNum = max(movNumsList);
end

if sum(renamedInds)<length(renamedInds)
    renameFlag = true;
else
    renameFlag = false;
end

%% rename cines if needed
if renameFlag

    datetimes = cell(1,length(camNamesList));
    cineDirs = cell(1,length(camNamesList));
    
    for camInd = 1:length(camNamesList)
        currCamInds = strcmp(camNamesList{camInd},{names_struct.camName});
        cineDirs{camInd} = names_struct(currCamInds);
    
            for cineFileInd = 1:sum(currCamInds)
                currName = cineDirs{camInd}(cineFileInd);
                datetime_str = [currName.year,currName.month,currName.day,...
                                currName.hour,currName.minute,currName.second(1:7)];
                cineDirs{camInd}(cineFileInd).triggerTime = datetime(datetime_str,...
                        'InputFormat','yyyyMMddHHmmss.SSSS');
            end
            datetimes{camInd} = [cineDirs{camInd}.triggerTime];
    end
    
    % compare datetimes and if a set of 3 vids is within the tol, group them
    for ind1 = 1:length(datetimes{1})
        [timeDiffs_2,ind2] = min(abs((datetimes{1}(ind1)-datetimes{2})));
        [timeDiffs_3,ind3] = min(abs((datetimes{1}(ind1)-datetimes{3})));
    
        tol = duration(0,0,0.1); % tolerance is 0.1 seconds (double triggering)
        if timeDiffs_2 < tol && timeDiffs_3 < tol
            movNum = movNum + 1;
            movNumStr = num2str(movNum,'%03.f');
    
            origFilenames = {cineDirs{1}(ind1).filename,cineDirs{2}(ind2).filename,...
                    cineDirs{3}(ind3).filename};
            origPaths = fullfile(pathToWatch,origFilenames);
            renameFilenames = strcat(camNamesList,'_',movNumStr,'.cine');
            cinePaths = fullfile(pathToWatch,renameFilenames);
    
            % rename files
            movefile(origPaths{1},cinePaths{1}); movefile(origPaths{2},cinePaths{2});
            movefile(origPaths{3},cinePaths{3});
    
            % also rename the .xml files
            [~,camNames,~] = fileparts(origPaths);
    
            xmlFilenames = strcat(camNamesList,'_',movNumStr,'.xml');
            xmlPaths = fullfile(pathToWatch,xmlFilenames);
    
            movefile(fullfile(pathToWatch,[camNames{1},'.xml']),xmlPaths{1});
            movefile(fullfile(pathToWatch,[camNames{2},'.xml']),xmlPaths{2});
            movefile(fullfile(pathToWatch,[camNames{3},'.xml']),xmlPaths{3});
    else
            continue
        end
    end
end

% get movie number list after renaming
cineDir = dir(fullfile(pathToWatch,'*.cine'));
renamed_out = regexp({cineDir.name},fnExpMovNum,'names');
renamedInds = ~(cellfun(@isempty,renamed_out));

if any(renamedInds)
    movsRenamed = renamed_out(renamedInds);
    renamedStruct = [movsRenamed{:}];
    movNumsList = str2double({renamedStruct.movieNum});
else
    disp('No movies to run reconstruction on')
    return
end

%% run analysis on renamed cines
for currMovNum = 8:12%unique(movNumsList)
    % check if there is triplet
    movNumStr = num2str(currMovNum,'%03.f');
    tripletCheck = sum(movNumsList==currMovNum);
    if tripletCheck ~= length(camNamesList)
        disp(['No matching triplet for movie ',movNumStr])
        continue
    else
        disp(['Found triplet for movie ',movNumStr,', running reconstruction'])
        analyzeOneFlyMovie(pathToWatch,currMovNum)
        disp(['Done analyzing movie ',movNumStr])
    end
end








