% Script to run analysis on renamed files. Need to modify the cnie2sparse
% code to make everything work.

pathToWatch = 'C:\Users\Abby\Cornell\Abby_test_tethered_code\01_23082024\' ; 
pathStruct = generatePathStruct(pathToWatch) ;
ExprNum = pathStruct.ExprNum;

clustFlag = true ; % shich version of analysis script to run
largePertFlag = false  ; % is it a large perturbation?
removeLegsFlag = false ; % try to remove legs in binary threshold?
alignBBoxFlag = false ; % try to align images to avoid clipping?
stopWingsFlag = false;

% 
% flyAnalysisMain(movNum, ExprNum, pathStruct, clustFlag, ...
%                 largePertFlag, removeLegsFlag) ;

% indexing for cameras
XZ = 2 ;
XY = 3 ;
YZ = 1 ;
movNum = 1;
movNumStr = num2str(movNum,'%03.f');
cineSuffix = ['_',movNumStr,'.cine'];
movieNum = movNumStr;

camNames = {'yz','xz','xy'};
camFilenames = strcat(camNames,cineSuffix);
cinFilenames = fullfile(pathToWatch,camFilenames);

twoFlies = 0;

load(fullfile(pathStruct.calibration,'calibration_easyWandData'))
DLT_matrix_CSV_filename = fullfile(pathStruct.calibration,...
    'calibration_dltCoefs.csv') ;

dlt_matrix = load(DLT_matrix_CSV_filename);

savePath = pathStruct.save;
if ~isfolder(savePath)
    mkdir(savePath)
end

prefixStr = ['mov_',movieNum];

movieFolder = fullfile(savePath,prefixStr);

if ~isfolder(movieFolder)
    mkdir(movieFolder)
end
hullFigPath = movieFolder;


%% FIND BACKGROUND ETC.
allBGcell = findBGTethered(pathToWatch);

% estimations for the time the fly comes in and out of the FOV of each
% camera. this part of the automation can be improved. check it or just set
% the "tin" and "tout" manually later.
% allTin  = zeros(3,1) ;
% allTout = zeros(3,1) ;
allTin = [-200;-200;-200];
allTout = [1400;1400;1400];

tin = -200;
tout = 1400;

% do this properly with the meta data later?
xcm = ones(2301,1);
ycm = ones(2301,1);

allXcm = {xcm*300,xcm*320,xcm*320};
allYcm = {ycm*340,ycm*350,ycm*400};
%---------------------------------------------
% save these results in case of error later?
% if savePointFlag
%    bg_savename = fullfile(savePath, prefixStr, 'BG.mat') ;
%    save(bg_savename,'allBGcell','allXcm','allYcm','allTin',...
%        'allTout','tin','tout')
% end
%  -----------------------------------------------------------------------
%% PERFORM BINARY THRESHOLDING ON IMAGES
%  -----------------------------------------------------------------------
% load phantom library
LoadPhantomLibraries();
RegisterPhantom(true);

cam = XY ;
try
    [all_fly_bw_xy, body_only_bw_xy, all_fly_thresholds_xy, xcm_xy, ycm_xy,...
        allAxlim_xy, DELTA, with_legs_bw_xy] = ...
        binaryThreshold(allBGcell{cam} , cinFilenames{cam}, tin,...
         tout, twoFlies, allXcm{cam}, allYcm{cam}, removeLegsFlag, ...
         stopWingsFlag) ;
catch exception
    msg = strcat('Error doing xy binary threshold for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    disp(msg)
    % fileID = fopen(errorPath,'a+') ;
    % fprintf(fileID, '%s\r\n', msg) ;
    % fclose(fileID) ;
    errorFlag = true ;
    return
end

cam = XZ ;
try
    [all_fly_bw_xz, body_only_bw_xz, all_fly_thresholds_xz, xcm_xz, ycm_xz,...
        allAxlim_xz, DELTA, with_legs_bw_xz] = ...
        binaryThreshold( allBGcell{cam}  , cinFilenames{cam}, tin,...
         tout, twoFlies, allXcm{cam}, allYcm{cam}, removeLegsFlag, ...
         stopWingsFlag) ;
catch exception
    msg = strcat('Error doing xz binary threshold for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    % fileID = fopen(errorPath,'a+') ;
    % fprintf(fileID, '%s\r\n', msg) ;
    % fclose(fileID) ;
    errorflag = true ;
    return
end

cam = YZ ;
try
    [all_fly_bw_yz, body_only_bw_yz, all_fly_thresholds_yz, xcm_yz, ycm_yz,...
        allAxlim_yz, DELTA, with_legs_bw_yz] = ...
        binaryThreshold( allBGcell{cam}  , cinFilenames{cam}, tin,...
         tout, twoFlies, allXcm{cam}, allYcm{cam}, removeLegsFlag, ...
         stopWingsFlag) ;
catch exception
    msg = strcat('Error doing yz binary threshold for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    disp(msg)
    % fileID = fopen(errorPath,'a+') ;
    % fprintf(fileID, '%s\r\n', msg) ;
    % fclose(fileID) ;
    errorFlag = true ;
    return
end
UnregisterPhantom();
UnloadPhantomLibraries();
%  -----------------------------------------------------------------------
%% COMBINE all_fly_bw_** INTO ONE STRUCTURE
%   (use frames DELTA+1 until Nimages-DELTA)
%  -----------------------------------------------------------------------

dim = max([all_fly_bw_xy.dim ; all_fly_bw_xz.dim ; all_fly_bw_yz.dim]) ; 
%dim = all_fly_bw_xy.dim ;
newdim = dim ;
% newdim(1) = dim(1) - 2*DELTA ;
% newdim(2) = 3 ; % three cams
newdim(2) = dim(2) - 2*DELTA ;
newdim(1) = 3 ; % three cams

all_fly_bw = init4D(newdim) ; 
% Nimages = newdim(1) ;
Nimages = newdim(2) ;
disp('Combining...')
for k=1:Nimages
    % ------------------
    % combine XY
    i1=XY ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize XY if need be
    bw_xy = getImage4D(all_fly_bw_xy, 1, i2+DELTA) ;
    if size(bw_xy,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_xy,1))/2) ;
        bw_xy = padarray(bw_xy,[pad_height,0],0,'both') ; 
    end
    if size(bw_xy,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_xy,2))/2) ;
        bw_xy = padarray(bw_xy,[0, pad_width],0,'both') ; 
    end
    all_fly_bw.mat(ind1vec, ind2vec) = bw_xy ;   
    
    % ------------------
    % combine XZ
    i1=XZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize XZ if need be
    bw_xz =  getImage4D(all_fly_bw_xz, 1, i2+DELTA) ; 
    if size(bw_xz,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_xz,1))/2) ;
        bw_xz = padarray(bw_xz,[pad_height,0],0,'both') ; 
    end
    if size(bw_xz,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_xz,2))/2) ;
        bw_xz = padarray(bw_xz,[0, pad_width],0,'both') ; 
    end
    all_fly_bw.mat(ind1vec, ind2vec) = bw_xz ;
    
    % ------------------
    % combine YZ
    i1=YZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize YZ if need be
    bw_yz =   getImage4D(all_fly_bw_yz, 1, i2+DELTA) ;   
    if size(bw_yz,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_yz,1))/2) ;
        bw_yz = padarray(bw_yz,[pad_height,0],0,'both') ; 
    end
    if size(bw_yz,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_yz,2))/2) ;
        bw_yz = padarray(bw_yz,[0, pad_width],0,'both') ; 
    end
    all_fly_bw.mat(ind1vec, ind2vec) = bw_yz ;
end

% ------------------------------------------------
% combine body_only_bw_** into one structure

dim = max([body_only_bw_xy.dim ; body_only_bw_xz.dim ; ...
    body_only_bw_yz.dim]) ; 
% dim = body_only_bw_xy.dim ;
newdim = dim ;
% newdim(1) = dim(1) - 2*DELTA ;
% newdim(2) = 3 ; % three cams
% body_only_bw = init4D(newdim) ; 
% Nimages = newdim(1) ;
newdim(2) = dim(2) - 2*DELTA ;
newdim(1) = 3 ; % three cams
body_only_bw = init4D(newdim) ; 
Nimages = newdim(2) ;

% WHEN dealing with body-only need to handle fucking delta.

for k=1:Nimages
    % --------------------
    % combine XY
    i1=XY ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize XY if need be
    bw_xy = getImage4D(body_only_bw_xy, 1, i2+DELTA) ;
    if size(bw_xy,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_xy,1))/2) ;
        bw_xy = padarray(bw_xy,[pad_height,0],0,'both') ; 
    end
    if size(bw_xy,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_xy,2))/2) ;
        bw_xy = padarray(bw_xy,[0, pad_width],0,'both') ; 
    end   
    body_only_bw.mat(ind1vec, ind2vec) = bw_xy ;   
    
    % --------------------
    % combine XZ
    i1=XZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize XZ if need be
    bw_xz = getImage4D(body_only_bw_xz, 1, i2+DELTA) ; 
    if size(bw_xz,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_xz,1))/2) ;
        bw_xz = padarray(bw_xz,[pad_height,0],0,'both') ; 
    end
    if size(bw_xz,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_xz,2))/2) ;
        bw_xz = padarray(bw_xz,[0, pad_width],0,'both') ; 
    end
    body_only_bw.mat(ind1vec, ind2vec) = bw_xz;   
    
    % --------------------
    % combine YZ
    i1=YZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize YZ if need be
    bw_yz = getImage4D(body_only_bw_yz, 1, i2+DELTA)  ;   
    if size(bw_yz,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_yz,1))/2) ;
        bw_yz = padarray(bw_yz,[pad_height,0],0,'both') ; 
    end
    if size(bw_yz,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_yz,2))/2) ;
        bw_yz = padarray(bw_yz,[0, pad_width],0,'both') ; 
    end
    body_only_bw.mat(ind1vec, ind2vec) = bw_yz ;   
end

% ------------------------------------------------------
% combine binarized images that include legs, if using

if isfield(with_legs_bw_xy, 'dim')
    dim = max([with_legs_bw_xy.dim ; with_legs_bw_xz.dim ; with_legs_bw_yz.dim]) ; 
    % dim = with_legs_bw_xy.dim ;
    newdim = dim ;
    % newdim(1) = dim(1) - 2*DELTA ;
    % newdim(2) = 3 ; % three cams
    newdim(2) = dim(2) - 2*DELTA ;
    newdim(1) = 3 ; % three cams
    
    with_legs_bw = init4D(newdim) ;
    % Nimages = newdim(1) ;
    Nimages = newdim(2) ;
    for k=1:Nimages
        % -------------------------
        % combine XY
        i1=XY ; i2=k ;
        ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
        ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
        % resize XY if need be
        bw_xy = getImage4D(with_legs_bw_xy, 1, i2+DELTA) ;
        if size(bw_xy,1) ~= dim(3)
            pad_height = round((dim(3) - size(bw_xy,1))/2) ;
            bw_xy = padarray(bw_xy,[pad_height,0],0,'both') ;
        end
        if size(bw_xy,2) ~= dim(4)
            pad_width = round((dim(4) - size(bw_xy,2))/2) ;
            bw_xy = padarray(bw_xy,[0, pad_width],0,'both') ;
        end
        with_legs_bw.mat(ind1vec, ind2vec) = bw_xy ;
        
        % -------------------------
        % combine XZ
        i1=XZ ; i2=k ;
        ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
        ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
        % resize XZ if need be
        bw_xz = getImage4D(with_legs_bw_xz, 1, i2+DELTA) ;
        if size(bw_xz,1) ~= dim(3)
            pad_height = round((dim(3) - size(bw_xz,1))/2) ;
            bw_xz = padarray(bw_xz,[pad_height,0],0,'both') ;
        end
        if size(bw_xz,2) ~= dim(4)
            pad_width = round((dim(4) - size(bw_xz,2))/2) ;
            bw_xz = padarray(bw_xz,[0, pad_width],0,'both') ;
        end
        with_legs_bw.mat(ind1vec, ind2vec) = bw_xz ;
        
        % -------------------------
        % combine YZ
        i1=YZ ; i2=k ;
        ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
        ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
        % resize YZ if need be
        bw_yz = getImage4D(with_legs_bw_yz, 1, i2+DELTA) ;
        if size(bw_yz,1) ~= dim(3)
            pad_height = round((dim(3) - size(bw_yz,1))/2) ;
            bw_yz = padarray(bw_yz,[pad_height,0],0,'both') ;
        end
        if size(bw_yz,2) ~= dim(4)
            pad_width = round((dim(4) - size(bw_yz,2))/2) ;
            bw_yz = padarray(bw_yz,[0, pad_width],0,'both') ;
        end
        with_legs_bw.mat(ind1vec, ind2vec) = bw_yz ;
    end
else
    with_legs_bw = [] ;
end

disp(['Done combining for movie ' movieNum])
%--------------------------------------------------------------------------
%% COMBINE CENTER-OF-MASS COORDINATES FOR EACH FRAME/CAMERA
%--------------------------------------------------------------------------
% CM_pos is the center-of-mass of the body in each image. 
% dimension of CM_pos is (3cameras, Nimages, 2coordinates)

% first adjust for any size changes we may have made to images
pad_amt_yz = all_fly_bw.dim(3:4) - all_fly_bw_yz.dim(3:4) ; 
pad_amt_xz = all_fly_bw.dim(3:4) - all_fly_bw_xz.dim(3:4) ; 
pad_amt_xy = all_fly_bw.dim(3:4) - all_fly_bw_xy.dim(3:4) ; 

% then combine CM measured from each camera into CM_pos
CM_pos = zeros(3, Nimages, 2) ;
ind = (1:Nimages) + DELTA ;
CM_pos(XY, :,1) = xcm_xy(ind) + pad_amt_xy(2)/2 ;
CM_pos(XY, :,2) = ycm_xy(ind) + pad_amt_xy(1)/2;
CM_pos(XZ, :,1) = xcm_xz(ind) + pad_amt_xz(2)/2;
CM_pos(XZ, :,2) = ycm_xz(ind) + pad_amt_xz(1)/2;
CM_pos(YZ, :,1) = xcm_yz(ind) + pad_amt_yz(2)/2;
CM_pos(YZ, :,2) = ycm_yz(ind) + pad_amt_yz(1)/2;


%% DEFINE PARAMS
% --------------

params.CAMERAS=[1 2 3];
params.NCAMS=3;
params.YZ = YZ ; 
params.XY = XY ;
params.XZ = XZ ;
params.fps = 8000 ;
% params.camerasPos=[easyWandData.DLTtranslationVector(:,:,2)';...
% easyWandData.DLTtranslationVector(:,:,1)';...
% easyWandData.DLTtranslationVector(:,:,3)'];  % row1=yz, row2=xz, row3=xy
params.cameraNames = ['yz';'xz';'xy'];
params.detectorLengthPix = all_fly_bw.dim(3:4) ;  % [imageHeight, imageWidth]
params.voxelSize = 50e-6 ; % 50 microns
%params.N=120; % 
params.volLength= 8e-3 ; % size of the square sub-vol cube to reconstruct (meters)
%params.voxelSize=0.4/120;
params.volCenter=[0,0,0];
params.focusPix=easyWandData.focalLengths;
%params.offsetsMatrix=[];
params.startTrackingTime   =  tin+DELTA  ; % plus 1 removed.
params.endTrackingTime     =  tout-DELTA ;
params.firstTrackableFrame =  tin+DELTA  ; % plus 1 removed.

% this is in fact "voxels per cm".
params.pixPerCM            = 350 ; % 232 ; % effective value. need to change that to be consistent with real units
% probably we should start from the wing legnth in cm and calculate how
% many voxels are in 1 wing length by dividing wingLcm/params.voxelSize

%---------------------------------------------
% save these results in case of error later?
% if savePointFlag
%    binaryThresh_savename = fullfile(savePath, prefixStr, 'binaryThresh.mat') ;
%    save(binaryThresh_savename, 'all_fly_bw_xy',  'body_only_bw_xy',...
%     'all_fly_thresholds_xy',  'xcm_xy', 'ycm_xy', 'allAxlim_xy', 'DELTA', ...
%     'all_fly_bw_xz', 'body_only_bw_xz', 'all_fly_thresholds_xz', 'xcm_xz',...
%     'ycm_xz', 'allAxlim_xz', 'all_fly_bw_yz', 'body_only_bw_yz', ...
%     'all_fly_thresholds_yz', 'xcm_yz', 'ycm_yz', 'allAxlim_yz', ...
%     'with_legs_bw_xy', 'with_legs_bw_xz', 'with_legs_bw_yz', 'all_fly_bw',...
%     'body_only_bw', 'with_legs_bw', 'params', 'easyWandData')
% end
%% VIZ body vs. whole fly featuring
% ---------------------------------
%{
figure('position',[ 94   584   560   160]);

ww = [-1 1 -1 1] * 48  ;

for k=1:Nimages
    for cam=1:3
        fly = getImage4D(all_fly_bw,cam,k) ;
        bod = getImage4D(body_only_bw, cam,k) ;
        rgb = zeros(512,512,3,'uint8') ;        
        rgb(:,:,1) = uint8(bod)*255 ; 
        rgb(:,:,2) = uint8(fly)*255 ; 
        subplot(1,3,cam) ; 
        imshow(rgb) ;
        xc = squeeze(CM_pos(cam,k,1)) ;
        yc = squeeze(CM_pos(cam,k,2)) ;
        axis([xc xc yc yc]+ww) ;
    end
    title(k) ;
    %saveas(gcf,['.\tmp\featuring_' num2str(k) '.png']) ;
    pause (0.05);
end
%}
%--------------------------------------------------------------------------
%% NEW 3D RECONSTRUCTION
%--------------------------------------------------------------------------
tic ;
disp('Doing hull reconstruction...')
try
    [ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
        wing1Res, wing1FrameStartInd, wing1FrameEndInd, ...
        wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
        hullReconstruction_mk6(params, CM_pos, all_fly_bw, body_only_bw,...
        dlt_matrix, easyWandData,[2,1,3]);
catch exception
    msg = strcat('Error reconstructing hulls for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    disp(msg)
    % fileID = fopen(errorPath,'a+') ;
    % fprintf(fileID, '%s\r\n', msg) ;
    % fprintf(fileID, '%s\r\n', ' ') ;
    % fclose(fileID) ;
    % errorFlag = true ;
    return
end
% [ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
%     wing1Res, wing1FrameStartInd, wingFrameEndInd, ...
%     wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
%     hullReconstruction_mk6(params, CM_pos, all_fly_bw, body_only_bw, dlt_matrix, easyWandData,[2,1,3]);
thull = toc ;

delete(gcp) ;

disp(['done calculating Hulls for movie ' movieNum]) ;

%---------------------------------------------
% save these results in case of error later?
% if savePointFlag
%    hullRecon_savename = fullfile(savePath, prefixStr, 'hullRecon.mat') ;
%    save(hullRecon_savename, 'bodyRes', 'bodyFrameStartInd', ...
%        'bodyFrameEndInd', 'wing1Res', 'wing1FrameStartInd',...
%        'wing1FrameEndInd', 'wing2Res', 'wing2FrameStartInd',...
%        'wing2FrameEndInd')
% end
%--------------------------------------------------------------------------
%% ANALYZE VOXEL RECONSTRUCTION
%--------------------------------------------------------------------------
t1 = clock ;

%diaryFile = ['myDiary.txt'] ;
% diaryFile = fullfile(savePath,prefixStr,'myDiary.txt');
% try
%     dos(['del ' diaryFile ]) ;
% catch
%     disp('Diary file does not exist.') ;
%     pause(2) ;
% end
% diary(diaryFile) ;
disp('Doing hull analysis...') 
plotHullFlag = false;
saveHullFigFlag = false;
try
    data = hullAnalysis_mk3 (bodyRes, wing1Res, wing2Res, params, ...
        mergedWingsFlag, [], 'test', plotHullFlag, saveHullFigFlag, ...
        hullFigPath);
catch exception
    msg = strcat('Error analyzing hulls for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    disp(msg)
    % fileID = fopen(errorPath,'a+') ;
    % fprintf(fileID, '%s\r\n', msg) ;
    % fprintf(fileID, '%s\r\n', ' ') ;
    % fclose(fileID) ;
    % errorFlag = true ;
    return
end
disp(['Done with hull analysis for movie ' movieNum])
t2 = clock ;
dt12 = t2 - t1 ;

% diary off ;
%--------------------------------------------------------------------------
%% SAVE RESULTS
%--------------------------------------------------------------------------


resultsFileName = [prefixStr,'_results'];

savePathFull = fullfile(savePath, prefixStr, [resultsFileName '.mat']) ; 
if exist('allBG','var')
save(savePathFull, 'data', 'bodyRes', 'bodyFrameStartInd', 'bodyFrameEndInd', ...
    'wing1Res', 'wing1FrameStartInd', 'wing1FrameEndInd', ...
    'wing2Res', 'wing2FrameStartInd', 'wing2FrameEndInd', 'mergedWingsFlag', ...
    'params',  'CM_pos', 'all_fly_bw', 'body_only_bw', 'with_legs_bw',...
    'dlt_matrix', 'easyWandData', ...
    'cinFilenames', 'allBG', 'Nimages', 'all_fly_bw_xy',  'body_only_bw_xy',...
    'all_fly_thresholds_xy',  'xcm_xy', 'ycm_xy', 'allAxlim_xy', 'DELTA', ...
    'all_fly_bw_xz', 'body_only_bw_xz', 'all_fly_thresholds_xz', 'xcm_xz',...
    'ycm_xz', 'allAxlim_xz', 'all_fly_bw_yz', 'body_only_bw_yz', ...
    'all_fly_thresholds_yz', 'xcm_yz', 'ycm_yz', 'allAxlim_yz', ...
    'tin', 'tout', 'resultsFileName')   
else
save(savePathFull, 'data', 'bodyRes', 'bodyFrameStartInd', 'bodyFrameEndInd', ...
    'wing1Res', 'wing1FrameStartInd', 'wing1FrameEndInd', ...
    'wing2Res', 'wing2FrameStartInd', 'wing2FrameEndInd', 'mergedWingsFlag', ...
    'params',  'CM_pos', 'all_fly_bw', 'body_only_bw', 'with_legs_bw',...
    'dlt_matrix', 'easyWandData', ...
    'cinFilenames', 'allBGcell', 'Nimages', 'all_fly_bw_xy',  'body_only_bw_xy',...
    'all_fly_thresholds_xy',  'xcm_xy', 'ycm_xy', 'allAxlim_xy', 'DELTA', ...
    'all_fly_bw_xz', 'body_only_bw_xz', 'all_fly_thresholds_xz', 'xcm_xz',...
    'ycm_xz', 'allAxlim_xz', 'all_fly_bw_yz', 'body_only_bw_yz', ...
    'all_fly_thresholds_yz', 'xcm_yz', 'ycm_yz', 'allAxlim_yz', ...
    'tin', 'tout', 'resultsFileName')
end
%dos(['move/y results_temp.mat ' resultsFileName '.mat']) ;

return