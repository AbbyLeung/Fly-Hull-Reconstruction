%--------------------------------------------------------------------------
% SIMULATE PITCH AND ROLL ERROR DUE TO EASYWAND FRAME MISALIGNMENT
%
% For a given calibration, computes the correction rotation R_corr that
% maps the easyWand frame to the true lab frame (xy camera -> (0,0,-1)).
%
% Then simulates a grid of "true" AHat vectors (varying yaw, pitch, and
% roll in the lab frame), transforms each into the easyWand frame using
% R_corr^-1, and computes the resulting pitch and roll errors.
%
% Output:
%   Figure 1 - heatmap of pitch error (easyWand pitch - true pitch) over
%              the yaw x pitch grid (independent of roll)
%   Figure 2 - heatmap of roll error (easyWand roll - true roll) over the
%              yaw x pitch grid for rho = 0
%   Figure 3 - roll error vs true roll for a representative hover point
%--------------------------------------------------------------------------

clear ; close all ;

%% --- USER INPUTS ---------------------------------------------------------
calibrationPath = 'Z:\Abby\_Opto_Mechanical_Data\08_29052025\calibration' ;
order           = [1, 2, 3] ; % YZ=col2, XZ=col1, XY=col3

% grid of true lab-frame yaw (psi) and pitch (beta) values to simulate
psi_range_deg  = -180 : 5 : 180 ;  % yaw  (deg)
beta_range_deg =  -90 : 5 :  90 ;  % pitch (deg), full range

% roll range for sweep
rho_range_deg  =  -90 : 5 :  90 ;  % roll (deg)

% representative hover point for the roll-vs-rho plot (Fig 3)
hover_psi_deg  =   0 ;
hover_beta_deg =  45 ;
%--------------------------------------------------------------------------

%% --- load calibration and compute R_corr --------------------------------
dlt = load(fullfile(calibrationPath, 'calibration_dltCoefs.csv')) ; 
raw = load(fullfile(calibrationPath, 'calibration_easyWandData.mat')) ;
if isfield(raw, 'easyWandData')
    easyWandData = raw.easyWandData ;
else
    fields       = fieldnames(raw) ;
    easyWandData = raw.(fields{1}) ;
end

pp  = easyWandData.ppts ;
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

xy_optical_axis = Rxy(1:3,1:3)' * [0;0;1] ;
xy_optical_axis = xy_optical_axis / norm(xy_optical_axis) ;

fprintf('XY camera optical axis in easyWand frame: [%.4f, %.4f, %.4f]\n', ...
    xy_optical_axis) ;
fprintf('Angle from (0,0,-1): %.3f deg\n', ...
    acosd(dot(xy_optical_axis, [0;0;-1]))) ;

% rotation that maps easyWand frame -> lab frame
target = [0;0;-1] ;
ax     = cross(xy_optical_axis, target) ;
axNorm = norm(ax) ;
if axNorm < 1e-10
    R_corr = eye(3) ;
    fprintf('No correction needed.\n') ;
else
    ax     = ax / axNorm ;
    ang    = acos(dot(xy_optical_axis, target)) ;
    K_skew = [0, -ax(3), ax(2) ; ax(3), 0, -ax(1) ; -ax(2), ax(1), 0] ;
    R_corr = eye(3) + sin(ang)*K_skew + (1-cos(ang))*(K_skew^2) ;
    fprintf('Correction rotation angle: %.3f deg\n', rad2deg(ang)) ;
end

% R_corr maps easyWand -> lab, so lab -> easyWand is R_corr^T
R_to_ew = R_corr' ;

%% --- simulate pitch error over yaw x pitch grid -------------------------
% pitch error does not depend on roll, so we only need the 2D (psi, beta) grid
Npsi  = length(psi_range_deg) ;
Nbeta = length(beta_range_deg) ;
Nrho  = length(rho_range_deg) ;

pitch_error_deg = zeros(Nbeta, Npsi) ;

for i = 1:Nbeta
    beta_true = beta_range_deg(i) * pi/180 ;
    for j = 1:Npsi
        psi_true = psi_range_deg(j) * pi/180 ;

        % true AHat in lab frame
        AHat_true = [cos(beta_true)*cos(psi_true) ;
                     cos(beta_true)*sin(psi_true) ;
                     sin(beta_true)] ;

        % what the easyWand frame gives (rotate true lab -> easyWand)
        AHat_ew = R_to_ew * AHat_true ;

        % pitch as computed by current code (asin of z-component in easyWand)
        beta_ew = asind(AHat_ew(3)) ;

        pitch_error_deg(i,j) = beta_ew - beta_range_deg(i) ;
    end
end

%% --- simulate roll error over yaw x pitch grid (at rho = 0) -------------
roll_error_rho0_deg = zeros(Nbeta, Npsi) ;

for i = 1:Nbeta
    beta_true = beta_range_deg(i) * pi/180 ;
    for j = 1:Npsi
        psi_true = psi_range_deg(j) * pi/180 ;

        AHat_true = [cos(beta_true)*cos(psi_true) ;
                     cos(beta_true)*sin(psi_true) ;
                     sin(beta_true)] ;

        % roll vector: body y-axis expressed in lab frame
        % eulerRotationMatrix(psi, beta, rho=0)' * [0;1;0]
        M_true    = eulerRotationMatrix(psi_true, beta_true, 0) ;
        rollVec_lab = M_true' * [0;1;0] ;

        % transform to easyWand frame
        AHat_ew    = R_to_ew * AHat_true ;
        rollVec_ew = R_to_ew * rollVec_lab ;

        % easyWand pitch & yaw -> rotM_YP as current code does it
        psi_ew_rad = atan2(AHat_ew(2), AHat_ew(1)) ;
        beta_ew_rad = asin(AHat_ew(3)) ;
        rotM_YP_ew = eulerRotationMatrix(psi_ew_rad, beta_ew_rad, 0) ;

        % apply to roll vector, zero x, re-normalise
        rotyb = rotM_YP_ew * rollVec_ew ;
        rotyb(1) = 0 ;
        if norm(rotyb) < 1e-10
            roll_error_rho0_deg(i,j) = NaN ;
            continue
        end
        rotyb = rotyb / norm(rotyb) ;
        rho_ew = acosd(rotyb(2)) * sign(rotyb(3)) ;

        roll_error_rho0_deg(i,j) = rho_ew - 0 ;  % true rho = 0
    end
end

%% --- simulate roll error vs true roll at hover point --------------------
roll_error_vs_rho = zeros(Nrho, 1) ;

psi_hover  = hover_psi_deg  * pi/180 ;
beta_hover = hover_beta_deg * pi/180 ;
AHat_hover = [cos(beta_hover)*cos(psi_hover) ;
              cos(beta_hover)*sin(psi_hover) ;
              sin(beta_hover)] ;
AHat_ew_hover    = R_to_ew * AHat_hover ;
psi_ew_hover     = atan2(AHat_ew_hover(2), AHat_ew_hover(1)) ;
beta_ew_hover    = asin(AHat_ew_hover(3)) ;
rotM_YP_ew_hover = eulerRotationMatrix(psi_ew_hover, beta_ew_hover, 0) ;

for kk = 1:Nrho
    rho_true = rho_range_deg(kk) * pi/180 ;

    M_true      = eulerRotationMatrix(psi_hover, beta_hover, rho_true) ;
    rollVec_lab = M_true' * [0;1;0] ;
    rollVec_ew  = R_to_ew * rollVec_lab ;

    rotyb = rotM_YP_ew_hover * rollVec_ew ;
    rotyb(1) = 0 ;
    if norm(rotyb) < 1e-10
        roll_error_vs_rho(kk) = NaN ;
        continue
    end
    rotyb = rotyb / norm(rotyb) ;
    rho_ew = acosd(rotyb(2)) * sign(rotyb(3)) ;

    roll_error_vs_rho(kk) = rho_ew - rho_range_deg(kk) ;
end

%% --- Figure 1: pitch error heatmap (unchanged) --------------------------
figure('Name', 'Pitch error due to easyWand frame misalignment', ...
    'Position', [50 550 750 500]) ;

imagesc(psi_range_deg, beta_range_deg, pitch_error_deg) ;
axis xy ;
colormap(redblue(256)) ;
cb = colorbar ;
cb.Label.String = 'Pitch error: \beta_{easyWand} - \beta_{true}  (deg)' ;
clim_val = max(abs(pitch_error_deg(:))) ;
clim([-clim_val, clim_val]) ;
xlabel('True yaw \psi (deg)') ;
ylabel('True pitch \beta (deg)') ;
title(sprintf(['Pitch error from easyWand frame misalignment\n' ...
    'XY camera %.2f deg from (0,0,-1)'], ...
    acosd(dot(xy_optical_axis, [0;0;-1])))) ;
yline(45, 'w--', 'LineWidth', 1.5) ;

%% --- Figure 2: roll error heatmap at rho = 0 ----------------------------
figure('Name', 'Roll error at rho=0 due to easyWand frame misalignment', ...
    'Position', [820 550 750 500]) ;

imagesc(psi_range_deg, beta_range_deg, roll_error_rho0_deg) ;
axis xy ;
colormap(redblue(256)) ;
cb = colorbar ;
cb.Label.String = 'Roll error: \rho_{easyWand} - \rho_{true}  (deg)  [at \rho_{true}=0]' ;
clim_val = max(abs(roll_error_rho0_deg(:)), [], 'omitnan') ;
if clim_val > 0
    clim([-clim_val, clim_val]) ;
end
xlabel('True yaw \psi (deg)') ;
ylabel('True pitch \beta (deg)') ;
title(sprintf(['Roll error from easyWand frame misalignment (\rho_{true}=0)\n' ...
    'XY camera %.2f deg from (0,0,-1)'], ...
    acosd(dot(xy_optical_axis, [0;0;-1])))) ;
yline(45, 'w--', 'LineWidth', 1.5) ;

%% --- Figure 3: roll error vs true roll at hover point -------------------
figure('Name', 'Roll error vs true roll at hover point', ...
    'Position', [50 50 600 400]) ;

plot(rho_range_deg, roll_error_vs_rho, 'k-', 'LineWidth', 1.5) ;
yline(0, '--', 'Color', [0.5 0.5 0.5]) ;
xlabel('True roll \rho (deg)') ;
ylabel('\Delta\rho (deg)  [easyWand \minus true]') ;
title(sprintf('Roll error vs true roll\n\\psi=%.0f°, \\beta=%.0f° (hover point)', ...
    hover_psi_deg, hover_beta_deg)) ;
grid on ;

fprintf('\nRoll error at rho=0, psi=0, beta=45: %.3f deg\n', ...
    roll_error_rho0_deg(beta_range_deg==45, psi_range_deg==0)) ;
fprintf('Max |roll error| at rho=0: %.3f deg\n', ...
    max(abs(roll_error_rho0_deg(:)), [], 'omitnan')) ;

%% --- helper: symmetric red-blue colormap --------------------------------
function cmap = redblue(n)
% red-white-blue colormap, symmetric around zero
half = floor(n/2) ;
r    = [linspace(0.8, 1, half), linspace(1, 1,   n-half)] ;
g    = [linspace(0,   1, half), linspace(1, 0,   n-half)] ;
b    = [linspace(1,   1, half), linspace(1, 0.2, n-half)] ;
cmap = [r', g', b'] ;
end
