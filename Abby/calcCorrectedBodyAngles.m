%--------------------------------------------------------------------------
% Recalculate body yaw, pitch, and roll after correcting for the easyWand
% coordinate frame not being aligned with the true lab frame (gravity along
% -z). Computes R_corr by forcing the XY camera optical axis to point along
% (0,0,-1), then rotates AHat and rollVectors before recomputing angles.
%
% INPUTS:
%   data            - processed data struct (must contain AHat, rollVectors,
%                     rhoTimes, params)
%   calibrationPath - path to folder containing calibration_dltCoefs.csv
%                     and calibration_easyWandData.mat
%   largePertFlag   - boolean (default: false). When true, uses
%                     calcPitchLargePert on the corrected AHat to avoid
%                     gimbal lock / 180-deg yaw flips near vertical.
%
% OUTPUTS:
%   psi_corr        - corrected body yaw angles [Nimages x 1], degrees
%   beta_corr       - corrected body pitch angles [Nimages x 1], degrees
%   rho_corr        - corrected body roll angles [Nimages x 1], degrees
%   R_corr          - 3x3 rotation matrix applied (easyWand -> lab frame)
%   psi_corr_smooth - low-pass smoothed yaw [Nimages x 1], degrees
%   beta_corr_smooth- low-pass smoothed pitch [Nimages x 1], degrees
%   rho_corr_smooth - low-pass smoothed roll [Nimages x 1], degrees
%--------------------------------------------------------------------------
function [psi_corr, beta_corr, rho_corr, R_corr, ...
          psi_corr_smooth, beta_corr_smooth, rho_corr_smooth] = ...
    calcCorrectedBodyAngles(data, calibrationPath, largePertFlag)
%--------------------------------------------------------------------------
%% inputs
if ~exist('largePertFlag','var') || isempty(largePertFlag)
    largePertFlag = false ;
end
% XY camera is always column 3 in the DLT CSV for this rig
XY_col = 3 ;
RAD2DEG = 180 / pi ;

%--------------------------------------------------------------------------
%% load calibration data
dlt      = load(fullfile(calibrationPath, 'calibration_dltCoefs.csv')) ;
raw      = load(fullfile(calibrationPath, 'calibration_easyWandData.mat')) ;
if isfield(raw, 'easyWandData')
    easyWandData = raw.easyWandData ;
else
    fields       = fieldnames(raw) ;
    easyWandData = raw.(fields{1}) ;
end

%--------------------------------------------------------------------------
%% compute R_corr from XY camera optical axis
f   = easyWandData.focalLengths' ;
pp  = easyWandData.ppts ;

% focal length and principal point for XY camera (always column 3)
f_xy  = f(XY_col) ;
pp_xy = pp(2*XY_col-1 : 2*XY_col) ;

Kxy   = [f_xy, 0, pp_xy(1) ; 0, f_xy, pp_xy(2) ; 0, 0, 1] ;
dltxy = dlt(:, XY_col) ;
Axy   = [dltxy(1:4)' ; dltxy(5:8)' ; dltxy(9:11)', 1] ;
Rxy   = Kxy \ Axy ;
Rxy   = Rxy / max([norm(Rxy(1:3,1)), norm(Rxy(1:3,2)), norm(Rxy(1:3,3))]) ;

xy_optical_axis = Rxy(1:3,1:3)' * [0;0;1] ;
xy_optical_axis = xy_optical_axis / norm(xy_optical_axis) ;

target = [0;0;-1] ;
ax     = cross(xy_optical_axis, target) ;
axNorm = norm(ax) ;

if axNorm < 1e-10
    R_corr = eye(3) ;
    fprintf('calcCorrectedBodyAngles: no correction needed (axis already aligned).\n') ;
else
    ax     = ax / axNorm ;
    ang    = acos(dot(xy_optical_axis, target)) ;
    K_skew = [0, -ax(3), ax(2) ; ax(3), 0, -ax(1) ; -ax(2), ax(1), 0] ;
    R_corr = eye(3) + sin(ang)*K_skew + (1-cos(ang))*(K_skew^2) ;
    fprintf('calcCorrectedBodyAngles: correction angle = %.3f deg\n', ...
        ang * RAD2DEG) ;
end

%--------------------------------------------------------------------------
%% rotate AHat into corrected lab frame
AHat_corr = (R_corr * data.AHat')' ;  % Nimages x 3
Nimages   = data.Nimages ;

%--------------------------------------------------------------------------
%% corrected yaw and pitch (+ rotM_YP for roll)
if largePertFlag
    % frame-to-frame wind-up: avoids 180-deg flip when body is near vertical
    [beta_corr, psi_corr, rotM_YP, largePertFlag] = ...
        calcPitchLargePert(AHat_corr, false) ;
    beta_corr = beta_corr(:) ;
    psi_corr  = psi_corr(:) ;
else
    psi_corr  = RAD2DEG * atan2(AHat_corr(:,2), AHat_corr(:,1)) ;
    beta_corr = RAD2DEG * asin(AHat_corr(:,3)) ;

    rotM_YP = zeros(3, 3, Nimages) ;
    for k = 1:Nimages
        rotM_YP(:,:,k) = eulerRotationMatrix( ...
            psi_corr(k)  / RAD2DEG, ...
            beta_corr(k) / RAD2DEG, 0) ;
    end
end

%--------------------------------------------------------------------------
%% corrected roll
if (~isfield(data,'rhoTimes')) || isempty(data.rhoTimes)
    warning('calcCorrectedBodyAngles: rhoTimes not found; roll set to zero.') ;
    rho_corr = zeros(Nimages, 1) ;
    return
end

% rotate rollVectors into corrected lab frame
rollVectors_corr = (R_corr * data.rollVectors')' ;

% build time vector matching calcAnglesRaw_Sam convention
fps = data.params.fps ;
if isfield(data, 'startAnalysisTimeMS')
    startTime = data.startAnalysisTimeMS * fps / 1000 ;
    endTime   = data.endAnalysisTimeMS   * fps / 1000 ;
else
    startTime = data.params.startTrackingTime ;
    endTime   = data.params.endTrackingTime ;
end
t = (startTime:endTime) / fps ;

[smoothed_rho, ~, ~, ~, ~] = calcBodyRoll(data.rhoTimes, rollVectors_corr, ...
    t, rotM_YP, data.params, largePertFlag) ;

rho_corr = smoothed_rho(:) ;

%--------------------------------------------------------------------------
%% smooth corrected angles (same filter params as smoothBodyAngles)
smoothingParams = setSmoothingParams() ;

% yaw: unwrap first (same logic as smoothBodyAngles)
init_window  = 20 ;
psi_smooth   = psi_corr ;
psi_init     = nanmedian(psi_smooth(1:min(init_window, length(psi_smooth)))) ;
psi_smooth   = psi_smooth - psi_init ;
psi_smooth   = unwrap(psi_smooth) ;
for i = 1:length(psi_smooth)
    while psi_smooth(i) < -180 ; psi_smooth(i) = psi_smooth(i) + 360 ; end
    while psi_smooth(i) >  180 ; psi_smooth(i) = psi_smooth(i) - 360 ; end
end
psi_smooth = psi_smooth + psi_init ;

if largePertFlag
    psi_smooth = (180/pi) * unwrap((pi/180) * psi_corr) ;
end

psi_corr_smooth  = filterEulerAngle(psi_smooth,  smoothingParams.yaw_filt_lvl) ;
beta_corr_smooth = filterEulerAngle(beta_corr,   smoothingParams.pitch_filt_lvl) ;
rho_corr_smooth  = filterEulerAngle(rho_corr,    smoothingParams.roll_filt_lvl) ;

end
