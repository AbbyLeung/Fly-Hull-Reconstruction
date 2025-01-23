function [data] = calcFlyAngles(inputData,largePertFlag, plotFlag)
%ESTIMATETETHEREDANGLSE Summary of this function goes here
%   Detailed explanation goes here
data = inputData;
[rhoTimes, rollVectors] = estimateRollVector(inputData,largePertFlag) ;
data.rhoTimes = rhoTimes ;
data.rollVectors = rollVectors ;

% below is from calcAnglesMain
defineConstantsScript

% stuff for manual correction
if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end

% if (isfield(data,'correctionTime'))
%     correctionTime=data.correctionTime;
% elseif(isfield(data,'manualCorrRangeMS'))
%     correctionTime=data.manualCorrRangeMS;
% else
%     correctionTime = [-10, 50] ;
% end

%% calculate raw angles
[anglesLabFrame, anglesBodyFrame, t, ~, ~, ~, smoothed_rho, rho_t, ...
    rho_samp, rotM_YP, rotM_roll, largePertFlag ] = ...
        calcAnglesRaw_Sam(data, plotFlag,largePertFlag) ;

%% unwrap and spline smooth wing angles
%------------------------------------------------
% right stroke angle
phiR = -anglesBodyFrame(:, PHIR) ;
ignoreIndR = unique([find(isnan(phiR))' ignoreFrames]) ;
%phiR = phiR + 360 ;

for i = 1:length(phiR)
    while phiR(i) < -90
        phiR(i) = phiR(i) + 360 ;
    end
    while phiR(i) > 270
        phiR(i) = phiR(i) - 360 ;
    end
end

if (~isempty(ignoreIndR))
    phiR(ignoreIndR) = NaN ;
end

%-----------------------------------------------
% left stroke angle
phiL = +anglesBodyFrame(:, PHIL) ;
ignoreIndL = unique([find(isnan(phiL))'  ignoreFrames])  ;

for i = 1:length(phiL)
    while phiL(i) < -90
        phiL(i) = phiL(i) + 360 ;
    end
    while phiL(i) > 270
        phiL(i) = phiL(i) - 360 ;
    end
end

if (~isempty(ignoreIndL))
    phiL(ignoreIndL) = NaN ;
end

% hampel filter to remove outliers
[~, hampelR] = hampel(phiR, 7,2) ;
[~, hampelL] = hampel(phiL, 7,2) ;

phiR(hampelR) = NaN ;
phiL(hampelL) = NaN ;

%--------------------------------------------------------------------------
%% get wing flip times
[fwdFlipTimesR, backFlipTimesR, fwdFlipIndR, backFlipIndR, fwdFlipPhiR,...
    backFlipPhiR, badIndicesR] = findWingFlipTimes_mk3 (t, phiR, plotFlag);
[fwdFlipTimesL, backFlipTimesL, fwdFlipIndL, backFlipIndL, fwdFlipPhiL,...
    backFlipPhiL, badIndicesL] = findWingFlipTimes_mk3 (t, phiL, plotFlag);

phiR(badIndicesR) = NaN ;
phiL(badIndicesL) = NaN ;

%--------------------------------------------------------------------------
%% store data in structure
anglesBodyFrame(:,PHIR) = -phiR ;
anglesBodyFrame(:,PHIL) = phiL ;

data.fwdFlipTimesR = fwdFlipTimesR ;
data.fwdFlipIndR = fwdFlipIndR ;
data.backFlipTimesR = backFlipTimesR ;
data.backFlipIndR = backFlipIndR ;
data.fwdFlipTimesL = fwdFlipTimesL ;
data.fwdFlipIndL = fwdFlipIndL ;
data.backFlipTimesL = backFlipTimesL ;
data.backFlipIndL = backFlipIndL ;
%data.params.pulseLengthMS = pulseLengthMS ;
data.anglesBodyFrame = anglesBodyFrame ;
data.anglesLabFrame = anglesLabFrame ;

% 4/23/20: want to include additional info: rotation matrices, 
% largePertFlag, ExprNum, MovNum, smoothed angles
data.rotM_YP = rotM_YP ; 
data.rotM_roll = rotM_roll ; 
data.largePertFlag = largePertFlag ; 

% data.ExprNum = exprNum ; 
% data.MovNum = movNum ; 

%--------------------------------------------------------------------------
%% find stroke amplitude
phiRfwdFlip = fwdFlipPhiR' ;  % fnval(sp_phiR_low, fwdFlipTimesR) ;
phiRbckFlip = backFlipPhiR' ; % fnval(sp_phiR_low, backFlipTimesR) ;
phiLfwdFlip = fwdFlipPhiL' ;  % fnval(sp_phiL_low, fwdFlipTimesL) ;
phiLbckFlip = backFlipPhiL' ; % fnval(sp_phiL_low, backFlipTimesL) ;

% right wing
phir_comb = [phiRfwdFlip, phiRbckFlip] ;
phir_t    = [fwdFlipTimesR, backFlipTimesR] ;
rmat = [phir_t', phir_comb'] ;
rmat = sortrows(rmat,1) ;

phir_amp = abs(diff(rmat(:,2))) ;
phir_amp_t = rmat(1:end-1,1) + diff(rmat(:,1))/2 ;
mid_stroke_r = rmat(1:end-1,2) + diff(rmat(:,2))/2 ;

% left wing
phil_comb = [phiLfwdFlip, phiLbckFlip] ;
phil_t    = [fwdFlipTimesL, backFlipTimesL] ;
lmat = [phil_t', phil_comb'] ;
lmat = sortrows(lmat,1) ;

phil_amp = abs(diff(lmat(:,2))) ;
phil_amp_t = lmat(1:end-1,1) + diff(lmat(:,1))/2 ;
mid_stroke_l = lmat(1:end-1,2) + diff(lmat(:,2))/2 ;

% add stroke amplitude and associated times to struct
data.phil_amp_t = phil_amp_t ;
data.phir_amp_t = phir_amp_t ;
data.phiL_amp = phil_amp ;
data.phiR_amp = phir_amp ;

% -------------------------------------------------------------------------
%% smooth angles (body and wing)
% --------------------------------------------
% smooth wing angles (lab and body frames)
[~, smoothAnglesMatR_Lab, ~, ~, ~ ] = smoothWingAngles(data, 'R','Lab') ;
[~, smoothAnglesMatL_Lab, ~, ~, ~ ] = smoothWingAngles(data, 'L','Lab') ;
[~, smoothAnglesMatR_Body, ~, ~, ~ ] = smoothWingAngles(data, 'R','Body') ;
[~, smoothAnglesMatL_Body, ~, ~, ~ ] = smoothWingAngles(data, 'L','Body') ;
% make sure phiR is negative in body frame
if (mode(sign(smoothAnglesMatR_Body(1,:))) > 0)
    smoothAnglesMatR_Body(1,:) = -1.*smoothAnglesMatR_Body(1,:) ; 
end

% --------------------------------------------------------------------
% smooth body angles (just in lab frame -- not defined in body frame)
[pitchSmooth, yawSmooth, rollSmooth] = smoothBodyAngles(data,largePertFlag) ;

% ------------------------------------
% create arrays for smoothed angles
% NB: need to take transpose for wing angle mats
anglesLabFrameSmooth = [yawSmooth, pitchSmooth, smoothAnglesMatR_Lab', ...
    smoothAnglesMatL_Lab', rollSmooth] ; 
anglesBodyFrameSmooth = zeros(data.Nimages, 8);
anglesBodyFrameSmooth(:,[PHIR, THETAR, ETAR, PHIL, THETAL, ETAL]) = ...
    [smoothAnglesMatR_Body', smoothAnglesMatL_Body'] ; 
    
% --------------------------------------
% add to data struct
data.anglesLabFrameSmooth = anglesLabFrameSmooth ; 
data.anglesBodyFrameSmooth = anglesBodyFrameSmooth ; 

end

