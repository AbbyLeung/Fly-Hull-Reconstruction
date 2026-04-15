%--------------------------------------------------------------------------
% COMPARE BODY PITCH: easyWand frame vs. corrected lab frame
%
% Computes body pitch two ways for a single data file:
%   (1) Old method: asin(AHat(3)) directly in easyWand coordinates
%   (2) Corrected method: rotate easyWand frame so xy camera optical axis
%       aligns with (0,0,-1), then compute asin(AHat_rotated(3))
%
% USAGE:
%   Set dataFile and calibrationPath below, or leave empty to use UI pickers.
%--------------------------------------------------------------------------

clear ; close all ;

%% --- USER INPUTS ---------------------------------------------------------
dataFile        = 'Z:\Abby\_Opto_Mechanical_Data\22_12022026\Analysis\Pitch Up\Expr_22_mov_044\Expr_22_mov_044_results.mat' ; % full path to a *_test.mat or *_results.mat
calibrationPath = 'Z:\Abby\_Opto_Mechanical_Data\22_12022026\calibration' ; % full path to calibration folder
order           = [2, 1, 3] ; % default: YZ=col2, XZ=col1, XY=col3 in dlt
%--------------------------------------------------------------------------

%% --- load data -----------------------------------------------------------
if isempty(dataFile)
    [fname, fpath] = uigetfile('*.mat', 'Select a data file (_test or _results)') ;
    if isequal(fname, 0), disp('Cancelled.') ; return ; end
    dataFile = fullfile(fpath, fname) ;
end
tmp = load(dataFile) ;
if isfield(tmp, 'data')
    data = tmp.data ;
else
    fields = fieldnames(tmp) ;
    data   = tmp.(fields{1}) ;
end

%% --- load calibration ---------------------------------------------------
if isempty(calibrationPath)
    calibrationPath = uigetdir(pwd, 'Select calibration folder') ;
    if isequal(calibrationPath, 0), disp('Cancelled.') ; return ; end
end

dlt          = load(fullfile(calibrationPath, 'calibration_dltCoefs.csv')) ;
raw          = load(fullfile(calibrationPath, 'calibration_easyWandData.mat')) ;
if isfield(raw, 'easyWandData')
    easyWandData = raw.easyWandData ;
else
    fields       = fieldnames(raw) ;
    easyWandData = raw.(fields{1}) ;
end

%% --- extract xy camera optical axis from DLT ----------------------------
pp  = easyWandData.ppts ; % principal points
f   = easyWandData.focalLengths' ;
f   = f(:, order) ;
pp  = pp(:, [2*order(1)-1, 2*order(1), ...
             2*order(2)-1, 2*order(2), ...
             2*order(3)-1, 2*order(3)]) ;

Kxy   = [f(3), 0, pp(5) ; 0, f(3), pp(6) ; 0, 0, 1] ;
dltxy = dlt(:, order(3)) ;
Axy   = [dltxy(1:4)' ; dltxy(5:8)' ; dltxy(9:11)', 1] ;
Rxy   = Kxy \ Axy ;
Rxy   = Rxy / max([norm(Rxy(1:3,1)), norm(Rxy(1:3,2)), norm(Rxy(1:3,3))]) ;

% optical axis of xy camera in easyWand world frame
xy_optical_axis = Rxy(1:3,1:3)' * [0;0;1] ;
xy_optical_axis = xy_optical_axis / norm(xy_optical_axis) ;

fprintf('XY camera optical axis in easyWand frame: [%.4f, %.4f, %.4f]\n', ...
    xy_optical_axis) ;
fprintf('Angle from (0,0,-1): %.3f deg\n', ...
    acosd(dot(xy_optical_axis, [0;0;-1]))) ;

%% --- compute correction rotation ----------------------------------------
% Find R_corr such that R_corr * xy_optical_axis = [0;0;-1]
target = [0;0;-1] ;
ax     = cross(xy_optical_axis, target) ;
axNorm = norm(ax) ;

if axNorm < 1e-10
    R_corr = eye(3) ;
    fprintf('No correction needed — xy axis already aligned with (0,0,-1).\n') ;
else
    ax  = ax / axNorm ;
    ang = acos(dot(xy_optical_axis, target)) ;
    K      = [0, -ax(3), ax(2) ; ax(3), 0, -ax(1) ; -ax(2), ax(1), 0] ;
    R_corr = eye(3) + sin(ang)*K + (1-cos(ang))*(K^2) ;
    fprintf('Correction rotation angle: %.3f deg\n', rad2deg(ang)) ;
end

% sanity check
fprintf('XY axis after correction: [%.4f, %.4f, %.4f]\n', R_corr * xy_optical_axis) ;

%% --- compute old and corrected pitch ------------------------------------
AHat          = data.AHat ;           % Nimages x 3, in easyWand frame
AHat_corr     = (R_corr * AHat')' ;  % Nimages x 3, in corrected lab frame

pitch_old_deg  = asind(AHat(:,3)) ;
pitch_corr_deg = asind(AHat_corr(:,3)) ;
diff_deg       = pitch_corr_deg - pitch_old_deg ;

fps = data.params.fps ;
t   = ((0:data.Nimages-1) + data.params.startTrackingTime) / fps ;

fprintf('\nPitch difference (corrected - old):\n') ;
fprintf('  Mean:       %.3f deg\n', mean(diff_deg, 'omitnan')) ;
fprintf('  Std:        %.3f deg\n', std(diff_deg, 'omitnan')) ;
fprintf('  Max |diff|: %.3f deg\n', max(abs(diff_deg), [], 'omitnan')) ;

%% --- plot ---------------------------------------------------------------
figure('Name', 'Pitch: easyWand frame vs corrected lab frame', ...
    'Position', [100 100 900 550]) ;

subplot(2,1,1) ;
hold on ;
plot(t, pitch_old_deg,  'b',  'LineWidth', 1.2, 'DisplayName', 'Old (easyWand frame)') ;
plot(t, pitch_corr_deg, 'r--','LineWidth', 1.2, 'DisplayName', 'Corrected (lab frame)') ;
hold off ;
ylabel('Body pitch \beta (deg)') ;
xlabel('Time (s)') ;
legend('Location', 'best') ;
title('Body pitch: easyWand frame vs corrected lab frame') ;
grid on ;

subplot(2,1,2) ;
plot(t, diff_deg, 'k', 'LineWidth', 1.2) ;
yline(0, '--', 'Color', [0.5 0.5 0.5]) ;
ylabel('\Delta\beta (deg)  [corrected \minus old]') ;
xlabel('Time (s)') ;
title(sprintf('Pitch difference  (mean = %.3f deg,  max |diff| = %.3f deg)', ...
    mean(diff_deg, 'omitnan'), max(abs(diff_deg), [], 'omitnan'))) ;
grid on ;
