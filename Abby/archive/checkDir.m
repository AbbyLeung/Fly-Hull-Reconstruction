function [currDir] = checkDir(pathToWatch,cineFileSize)
%CHECKDIR Summary of this function goes here
%   Detailed explanation goes here
% check input(s)
if ~exist('pathToWatch','var') || isempty(pathToWatch)
    pathToWatch = 'I:\Fly_Experiment\antenna\001_17102024';
end

drivenames = {'I:','J:','K:','L:'};
pathSplit = split(pathToWatch,'\');
pathsToWatch = fullfile(drivenames,pathSplit{2:end});


% get search expression for video file names
% two digits in day
fnExp1 = ['(?<camName>Camera[1234])_Y(?<year>\d{4})',...
    '(?<month>\d{2})(?<day>\d{2})H(?<hour>\d{2})(?<minute>\d{2})', ...
    '(?<second>\d+.\d+)']; 
% renamed file format, such as xy_001
fnExp2 = '(?<camName>Camera[1234])_(?<movieNum>\d{3})';
% one digit in day AND ALSO a space between the month and day after PCC
% update
fnExp3 = ['(?<camName>Camera[1234])_Y(?<year>\d{4})',...
    '(?<month>\d{2}) (?<day>\d{1})H(?<hour>\d{2})(?<minute>\d{2})', ...
    '(?<second>\d+.\d+)']; 
fnExp = [fnExp1,'|',fnExp3]; 

% initialize global params for each camera
cam1EventCount = 0 ; cam2EventCount = 0 ; cam3EventCount = 0 ; cam4EventCount = 0;
cam1EventTimes = [] ; cam2EventTimes = [] ; cam3EventTimes = [] ; cam4EventTimes = [] ;
cam1FileNames = {} ;  cam2FileNames = {} ;  cam3FileNames = {} ; cam4FileNames = {} ;
CamTimes_list = {};
queue = [] ; %queue for video numbers that need analyzing
movFileExt = [] ;

% estimate cine file size to determine params for checking cine save status
cinFileSizeDefault = 1.0e9 ; 
if ~exist('cinFileSize','var') || isempty(cineFileSize)
    cineFileSize = cinFileSizeDefault ;
end

% Define cluster
try
    clust = parcluster('myLocalCluster') ;
catch
    clust = parcluster('local') ;
end
clust.JobStorageLocation = pathToWatch ;
clust.HasSharedFilesystem = true ;
analysis_job_cell = cell(1) ;
%mp4_job_cell = cell(1) ;

% get save path info
pathStruct = generatePathStruct(pathToWatch) ;
ExprNum = pathStruct.ExprNum ;
cd(pathToWatch)

% define the watcher for camera 1
fileObjXY = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjXY.EnableRaisingEvents = true ;
fileObjXY.Filter = 'Camera1*.cine' ;  % ['xy*' movFileExt] ;
addlistener(fileObjXY, 'Created', @(src, evt) onChange(src, evt)) ;

fprintf(['\n ' repmat('%',1,100) '\n' ])
fprintf('\n Running experiment: \n %s \n', pathToWatch)
fprintf(['\n ' repmat('%',1,100) '\n'])

totalMovies = 0;

while true
    % brief pause
    pause(1)

    % -----------------------------------------------------------------
    % if each file listener has seen a file be created, add this
    % triplet to queue
    if (cam1EventCount > 0) && (cam2EventCount > 0) && (cam3EventCount > 0) ...
            && (cam4EventCount)
        % give a little time (after files are created) to allow all movies
        % to be saved
        pause(10)
       
        % get filenames for matching triplet of movies (also check that
        % there are four movies with matching times)    
        [proceedFlag, cam1_fn, cam2_fn, cam3_fn, cam4_fn] = groupMovies(cam1EventTimes,...
            cam2EventTimes, cam3EventTimes, cam4EventTimes, ...
            cam1FileNames, cam2FileNames, cam3FileNames, cam4FileNames) ;
        %adds to queue once 4 videos (4 cameras) are captured
        if proceedFlag
            totalMovies = totalMovies+1;
            add2queue(cam1_fn, cam2_fn, cam3_fn, cam4_fn,totalMovies) ;
        else
            % do nothing
            fprintf('No matching movies to add to queue \n')
        end
    end

    % --------------------------------------------------------------
    % if we have something in the queue, run analysis
    if numel(queue) ~= 0 %do nothing if queue is empty
        % if we have something in queue, run analysis
        % make folder for movie
        analysis_job_cell{queue(1)+1} = batch(clust,'runFullAnalysis',...
            0,{exprNum,queue(1), pathsToWatch},'CaptureDiary',true,...
            'Pool', 1) ;

        % tag analysis job with movie number
        analysis_job_cell{queue(1)+1}.Tag = num2str(queue(1),'%03d') ;
       
        % after running analysis on most recent queue item (sent to batch),
        % remove from queue
        queue = queue(2:numel(queue)) ; %remove first element
        continue
    end
   
    % ----------------------------------------------------------------
    % delete completed jobs to free up cluster resources
    completedJobs = findJob(clust,'State','finished') ;
    if ~isempty(completedJobs)
        fprintf('[DONE] Completed analysis for movie %s (%s)\n', ...
            completedJobs.Tag, datetime("now"))
        delete(completedJobs) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALLBACK FUNCTIONS FOR FILE LISTENERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function onChange(~, evt)
    % -----------------------------------------------------
    % gets input from current file being saved, increments
    % movie-specific counters, and writes to log
    % -----------------------------------------------------
    % give time for a new file to show up in all 4 directories
    pause(15)
    % search filename for movie info
    eventName = char(evt.Name) ;
    out_names = regexp(eventName, fnExp,'names') ;
   
    % if we haven't already specified cine file extension, get it here
    if ~exist('movFileExt','var') || isempty(movFileExt)
        [~, ~, movFileExt] = fileparts(eventName) ;
    end
    % read out camera name (should be present regardless of file
    % format)
   
    camName = out_names.camName ;
    disp(['detected change in ',camName])
   

    % -------------------------------------------------------
    % get datenum information
    if ~isempty(out_names)
        trigger_datestr = strjoin({out_names.year, out_names.month, ...
            out_names.day, out_names.hour, out_names.minute,...
            out_names.second(1:6)},' ') ;
        trigger_datenum = datenum(trigger_datestr, ...
            'yyyy mm dd HH MM SS.FFF') ;
    else
        % otherwise, use current datetime as a proxy for trigger time
        % NB: THIS IS LAZY -- SHOULD GET FROM XML FILE
        trigger_datenum = datenum(datetime('now')) ;
    end

    % ----------------------------------
    % update log file
    change = strcat(camName, ': ', datestr(trigger_datenum)) ;
    fileID = fopen(strcat(pathToWatch,'\',eventName(1:2),'log.txt'),'a+') ;
    fprintf(fileID, '%s\r\n', change) ;
    fclose(fileID) ;

    % ---------------------------------
    % update global counters
    switch camName
        case 'Camera1'
            cam1FileNames{end+1} = eventName ;
            cam1EventTimes = [cam1EventTimes; trigger_datenum] ;
            cam1EventCount = cam1EventCount + 1 ;
        case 'Camera2'
            cam2FileNames{end+1} = eventName ;
            cam2EventTimes = [cam2EventTimes; trigger_datenum] ;
            cam2EventCount = cam2EventCount + 1 ;
        case 'Camera3'
            cam3FileNames{end+1} = eventName ;
            cam3EventTimes = [cam3EventTimes; trigger_datenum] ;
            cam3EventCount = cam3EventCount + 1 ;
        case 'Camera4'
            cam4FileNames{end+1} = eventName ;
            cam4EventTimes = [cam4EventTimes; trigger_datenum] ;
            cam4EventCount = cam4EventCount + 1 ;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
% take a set of three movies recently saved, make sure their
% timing matches up, and then output filenames for these movies
%
% NB: taking event times and filenames as inputs to reduce risk of
% list being updated during computation
% ---------------------------------------------------------------------
    function [proceedFlag, cam1_fn, cam2_fn, cam3_fn,cam4_fn] = groupMovies(...
            cam1Times,cam2Times, cam3Times, cam4Times,...
        cam1Names, cam2Names, cam3Names,cam4Names)
   
    % define datenum tolerance for matching triplet
    tol = 5e-5 ;
    % *** WAS 5e-4 (in days), currently ~4.32s

    % find all combinations of possible datenums in pre-queue
    camTimes_list = {cam1Times, cam2Times, cam3Times,cam4Times};
    [camTimes_list{end:-1:1}] = ndgrid(camTimes_list{end:-1:1});
    
    N_cams = numel(camTimes_list);
    eventTimeCombs = reshape(cat(N_cams,camTimes_list{:}),[],N_cams);

    % get combination with smallest difference in event timing
    timeDist = zeros(1,height(eventTimeCombs));

    for row = 1:height(eventTimeCombs)
        currComb = eventTimeCombs(row,:);
        timeDist(row) = max(abs(diff(nchoosek(currComb,2),[],2)));
    end
        
    [minTimeDist, minInd] = min(timeDist) ;

    % make sure the total distance between time points is within
    % tolerance
    if (minTimeDist < tol)
        % in this case, distance is within tolerance and we should
        % select this combination
        proceedFlag = true ;

        % find the indices in camera event times array corresponding towhosv
        % the datenum combination with minimal time gap
        [~, cam1_ind] = min(abs(cam1Times - eventTimeCombs(minInd,1))) ;
        [~, cam2_ind] = min(abs(cam2Times - eventTimeCombs(minInd,2))) ;
        [~, cam3_ind] = min(abs(cam3Times - eventTimeCombs(minInd,3))) ;
        [~, cam4_ind] = min(abs(cam4Times - eventTimeCombs(minInd,4))) ;

        % use these indices to get the matching filenames
        cam1_fn = cam1Names{cam1_ind} ;
        cam2_fn = cam2Names{cam2_ind} ;
        cam3_fn = cam3Names{cam3_ind} ;
        cam4_fn = cam4Names{cam4_ind} ;
    else
        proceedFlag = false ;
        cam1_fn = [] ; cam2_fn = [] ; cam3_fn = [] ; cam4_fn = [];
        decrementCounters()
    end
end

% ------------------------------------------------------------
% take 4 movie files and add them to the queue for
% analysis
%
% To add: merge logs
% -------------------------------------------------------------
    function add2queue(cam1_fn, cam2_fn, cam3_fn,cam4_fn, totalMovies)
    % -------------------------------------------------------------
    % make sure cin files are completely saved before adding to
    % queue (i.e. make sure file size is sufficiently large)
    pause(1) ;
    cam1Bytes = 0; cam2Bytes = 0; cam3Bytes = 0 ; cam4Bytes = 0;
    byteCounter = 0 ;
    while (cam1Bytes < cinFileSize || ...
            cam2Bytes < cinFileSize || ...
            cam3Bytes < cinFileSize || ...
            cam4Bytes < cinFileSize) && (byteCounter < byteCountMax)
        % give a little time to allow save to progress
        pause(5)
        % change the correct filenames back to wrong ones so that the
        % script is happy
        % get current directory for camera movie files
        cam1dir = dir(fullfile(pathsToWatch{1}, ['Camera1*' movFileExt])) ;
        cam2dir = dir(fullfile(pathsToWatch{2}, ['Camera2*' movFileExt])) ;
        cam3dir = dir(fullfile(pathsToWatch{3}, ['Camera3*' movFileExt])) ;
        cam4dir = dir(fullfile(pathsToWatch{4}, ['Camera4*' movFileExt])) ;


        % get index in directory for current movie triplet filenames
        cam1_ind = arrayfun(@(x) strcmp(x.name, cam1_fn), cam1dir) ;
        cam2_ind = arrayfun(@(x) strcmp(x.name, cam2_fn), cam2dir) ;
        cam3_ind = arrayfun(@(x) strcmp(x.name, cam3_fn), cam3dir) ;
        cam4_ind = arrayfun(@(x) strcmp(x.name, cam4_fn), cam4dir) ;

        % check current file size
        cam1Bytes = cam1dir(cam1_ind).bytes ;
        cam2Bytes = cam2dir(cam2_ind).bytes ;
        cam3Bytes = cam3dir(cam3_ind).bytes ;
        cam4Bytes = cam4dir(cam4_ind).bytes ;
       
        % increment counter to determine how many loops we've gone through
        byteCounter = byteCounter + 1 ;
    end
   
    % if we broke while loop because of too many loops, it might mean file
    % size doesn't match -- just skip over these for now and deal with them
    % later
    if byteCounter == byteCountMax
        divertData({cam1_fn, cam2_fn, cam3_fn,cam4_fn}, {movFileExt, '.xml'})
        return
    end
   
    % ... otherwise pause to allow any final saving processes to take place
    pause(5)

    % -----------------------------------------------------------------
    % if applicable, rename movies from datetime filename format to 3
    % digit integer format
    regexp_check = regexp({cam1_fn, cam2_fn, cam3_fn,cam4_fn}, fnExp, 'names') ;

    % first check if (unaltered) files have associated movie numbers
    % already. NB: numberMismatchFlag should correspond to cases wherein
    % the movie numbers don't match up; assignNumberFlag should correspond
    % to cases wherein the files are being saved with the trigger time in
    % the filename
    movNums = cellfun(@(y) str2double(y.movieNum), regexp_check, 'ErrorHandler',@movNumErrorFunc) ;
    %numberMismatchFlag = 0; %any(isnan(movNums)) | ...
    %    ~all(movNums(:) == movNums(1)) ;
    assignNumberFlag = true; % all(isnan(movNums)) ; %Made this always true, hopefully that doesn't cause a bug
   
    % regardless of reason why, define flag for any renaming that needs to
    % occur
    % renameFlag = assignNumberFlag | numberMismatchFlag ;
    renameFlag = true;
   
    % generate lists that we can loop over to make process of
    % renaming easier (if we're renaming)
    if renameFlag
        old_fn_list = {cam1_fn, cam2_fn, cam3_fn,cam4_fn} ;
        cam_names = {'Camera1', 'Camera2', 'Camera3','Camera4'} ;
        file_type_list = {movFileExt, '.xml'} ;
    end
   
    % if either we have non-matching movie numbers OR datetime format
    % movie filenames, rename the cine and xml files
    if assignNumberFlag
        % if we need to rename movies, want to get current maximum
        % movie number (to avoid overwriting)
        movieCountCurr = getMovieCounter() ; % current max movie number
       
        % TESTING -- with datetime file format, should have less risk of
        % overwriting, so just use (movieCountCurr + 1) as new movie number
        if totalMovies <= movieCountCurr
            movieCounter = movieCountCurr + 1 ;
        else
            movieCounter = totalMovies+1;
        end
      
    else
        % in this case, all movie numbers match, so we're good
        movieCounter = movNums(1) ;
    end
   
    % -------------------------------------------------------------
    % now that we have movie number sorted out (hopefully), perform
    % renaming if necessary
    if renameFlag
        for mm = 1:length(cam_names)
            % get current camera name
            cam = cam_names{mm} ;
            % get basename for old filename (for this camera)
            [~, base_name, ~] = fileparts(old_fn_list{mm}) ;
           
            % loop over file types
            for nn = 1:length(file_type_list)
                % current file extension
                ext_curr = file_type_list{nn} ;
               
                % rename file
                old_fn = fullfile(pathsToWatch{mm}, [base_name ext_curr]) ;
                new_fn = fullfile(pathsToWatch{mm}, ...
                    [cam '_' num2str(movieCounter,'%03d') ext_curr]) ;
               
                if ~strcmp(old_fn, new_fn)
                    status = movefile(old_fn, new_fn) ;
                    if ~status
                        fprintf('Error changing filename from %s to %s \n',...
                            old_fn, new_fn)
                        decrementCounters(cam1_fn, cam2_fn, cam3_fn, cam4_fn)
                        return
                    end
                end
            end
        end
    end

    % ----------------------------------------------------------------
    % add current movie to the queue and decrease event counters/remove
    % entries from event time lists
    queue(end+1) = movieCounter ;
    decrementCounters(cam1_fn, cam2_fn, cam3_fn, cam4_fn)
    fprintf('[ADD] Movie %03d added to analysis queue (%s) \n', ...
        movieCounter, datetime("now"))

end

% --------------------------------------------------------------------
% get maximum number of current movie in directory and return that
% number (or zero, if no movies)
function movCount = getMovieCounter()
    datadir = dir(fullfile(pathToWatch, ['*' movFileExt])) ;
    movNames = {datadir(:).name} ;
    if isempty(movNames)
        movCount = 0 ;
    else
%         searchExp = '(?<camName>[xyz]+)_(?<movieNum>\d+)' ;
        out_names = regexp(movNames, fnExp, 'names') ;
        movNums = cellfun(@(x) str2double(x.movieNum), out_names, 'ErrorHandler', @movNumErrorFunc) ;
        movCount = max([movNums, 0],[],'omitnan') ;
       
    end
end

% ---------------------------------------------------------------------
%% decrease event counters for the three cameras
    function decrementCounters(cam1_fn, cam2_fn, cam3_fn, cam4_fn)
    % decrease event counter by 1 for each movie
    cam1EventCount = cam1EventCount - 1 ;
    cam2EventCount = cam2EventCount - 1 ;
    cam3EventCount = cam3EventCount - 1 ;
    cam4EventCount = cam4EventCount - 1 ;
   
    % if we have specific files as inputs, remove the event times and file
    % names
    if (nargin > 0)
        % get indices for filenames
        cam1_ind = cellfun(@(x) strcmp(x, cam1_fn), cam1FileNames) ;
        cam2_ind = cellfun(@(x) strcmp(x, cam2_fn), cam2FileNames) ;
        cam3_ind = cellfun(@(x) strcmp(x, cam3_fn), cam3FileNames) ;
        cam4_ind = cellfun(@(x) strcmp(x, cam4_fn), cam4FileNames) ;
       
        % remove filenames and event numbers, decrement event counter
        cam1FileNames = cam1FileNames(~cam1_ind) ;  % xy
        cam1EventTimes = cam1EventTimes(~cam1_ind) ;
       
        cam2FileNames = cam2FileNames(~cam2_ind) ;  % xz
        cam2EventTimes = cam2EventTimes(~cam2_ind) ;
       
        cam3FileNames = cam3FileNames(~cam3_ind) ;  % yz
        cam3EventTimes = cam3EventTimes(~cam3_ind) ;  

        cam4FileNames = cam4FileNames(~cam4_ind) ; 
        cam4EventTimes = cam4EventTimes(~cam4_ind) ; 
    end
end

%% remove data from main directory due to potential save error
function divertData(fn_list, file_type_list)
    t_now = datetime('now') ;
    tempSaveDir = fullfile(pathToWatch, string(t_now,"yyyyMMddhhmmss")) ;
    mkdir(tempSaveDir)
   
    % loop over cameras and file types to move data
    for mm = 1:length(fn_list)
        % get basename for old filename (for this camera)
        [~, base_name, ~] = fileparts(fn_list{mm}) ;
        for nn = 1:length(file_type_list)
            % current file extension
            ext_curr = file_type_list{nn} ;
            fn_wExt = [base_name ext_curr] ;
           
            % move file
            old_fn = fullfile(pathToWatch, fn_wExt) ;
            new_fn = fullfile(tempSaveDir, fn_wExt) ;
            movefile(old_fn, new_fn) ;
        end
    end
   
    % after moving files, exit without adding to queue
    decrementCounters(fn_list{1}, fn_list{2}, fn_list{3}, fn_list{4})
end

%Function to handle errors thrown when a movie has no movie number
function movNum = movNumErrorFunc(S,varargin)
    %warning(S.identifier, S.message); 
    movNum = NaN;
end
end

