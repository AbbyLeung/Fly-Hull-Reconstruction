function hullInspectionGUI(resultsPath)
%HULLINSPECTIONGUI  Rapidly browse hull projections and flag bad frames.
%
%   hullInspectionGUI()             — opens file browser
%   hullInspectionGUI(resultsPath)  — loads file directly
%
%   Workflow
%   --------
%   1. Load a *_results.mat file (same folder as your cleaned data).
%   2. Use arrow keys (or buttons) to step through frames.
%   3. Press SPACE (or the Flag button) to mark a bad frame.
%   4. Press S to save flags → writes *_badFrames.mat alongside the data.
%   5. Open the correction GUI — it auto-loads the bad frames list so that
%      "Next Bad Frame" navigates only to flagged frames.
%
%   Keyboard shortcuts
%   ------------------
%   ←  / →        previous / next frame
%   Shift+← / →   jump back / forward 10 frames
%   ↑  / ↓        previous / next FLAGGED frame
%   Space         toggle flag on current frame
%   S             save flags
%   V             toggle vector overlay on/off
%
%   Voxel colour key
%   ----------------
%   White  = body
%   Red    = right wing
%   Blue   = left wing
%   Dim grey silhouette = original camera binary mask
%
%   Vector colour key (when a cleaned/corrected file is loaded)
%   ----------------------------------------------------------
%   Yellow  = body long axis (A-hat)
%   Red     = right-wing span      Orange = right-wing chord
%   Blue    = left-wing span       Cyan   = left-wing chord
%   Time display: t = (frame - 1 + params.startTrackingTime) / params.fps
%
%   Requires getImage4D and dlt_inverse to be on the MATLAB path
%   (both live in this project under core functions / camera calibration).
%
%   Camera / DLT mapping assumed:  order = [2 1 3]
%     YZ camera  → dlt_matrix(:,2)
%     XZ camera  → dlt_matrix(:,1)
%     XY camera  → dlt_matrix(:,3)
%   If your data used a different order, edit the s.dltCols line below.

if nargin < 1, resultsPath = ''; end

%% ── Figure ───────────────────────────────────────────────────────────────
FIGW = 1520;  FIGH = 840;
fig = figure( ...
    'Name',            'Hull Inspection GUI', ...
    'NumberTitle',     'off', ...
    'MenuBar',         'none', ...
    'ToolBar',         'none', ...
    'Color',           [0.94 0.94 0.94], ...
    'Position',        [20 40 FIGW FIGH], ...
    'KeyPressFcn',     @keyPressed, ...
    'CloseRequestFcn', @(src,~) delete(src));

%% ── State ────────────────────────────────────────────────────────────────
s.loaded        = false;
s.Nimages       = 0;
s.frameCurr     = 1;
s.flagged       = false(0,1);
s.resultsPath   = '';
s.data          = [];
s.dlt_matrix    = [];
s.all_fly_bw    = [];
s.body_only_bw  = [];
s.params        = [];
s.CM_pos        = [];
s.startIdx      = [];   % precomputed first row in data.res for each frame
s.endIdx        = [];   % precomputed last  row in data.res for each frame
s.imH           = 512;
s.imW           = 512;
% Camera index → DLT column mapping  (order = [2 1 3] convention)
% Row k: [ camera_index_in_getImage4D,  dlt_column ]
% YZ=cam1→col2,  XZ=cam2→col1,  XY=cam3→col3
s.camIdx        = [1  2  3];   % filled from params.YZ/XZ/XY on load
s.dltCols       = [2  1  3];   % dlt_matrix column for YZ / XZ / XY
s.camLabels     = {'YZ', 'XZ', 'XY'};
s.zoom          = 90;          % half-width of zoom window (pixels)

% Vector overlay (from cleaned/corrected data file)
s.showVectors    = true;       % toggle with V or button (body axis + span)
s.showChord      = false;      % toggle with C or button (chord vectors)
s.hasVectors     = false;      % populated in doLoad
s.bodyCM         = [];
s.rightWingCM    = [];
s.leftWingCM     = [];
s.AHat           = [];         % body long axis  (Nx3, unit vec in metres)
s.rightSpanHats  = [];
s.leftSpanHats   = [];
s.rightChordHats = [];
s.leftChordHats  = [];
s.vecScaleBody   = 0.0018;     % metres — visible length of body axis vector
s.vecScaleSpan   = 0.0022;     % metres — wing span vector
s.vecScaleChord  = 0.0008;     % metres — wing chord vector

% Time-axis
s.fps               = NaN;
s.startTrackingTime = NaN;     % absolute frame number that corresponds to t=0

%% ── Top bar ──────────────────────────────────────────────────────────────
BG  = [0.94 0.94 0.94];
BTN = [0.85 0.85 0.85];
FG  = [0.10 0.10 0.10];

uicontrol('Style','pushbutton','String','Load Results File', ...
    'Units','pixels','Position',[8 FIGH-42 160 30], ...
    'BackgroundColor',[0.22 0.48 0.22],'ForegroundColor','w', ...
    'FontSize',9,'FontWeight','bold','Callback',@loadFile);

s.hFileLabel = uicontrol('Style','text','String','No file loaded.', ...
    'Units','pixels','Position',[178 FIGH-40 FIGW-190 24], ...
    'HorizontalAlignment','left','FontSize',8, ...
    'BackgroundColor',BG,'ForegroundColor',[0.3 0.3 0.3]);

%% ── Three camera axes ────────────────────────────────────────────────────
axY = 0.185;  axH = 0.635;  axW = 0.305;
s.ax(1) = axes('Units','normalized','Position',[0.010 axY axW axH], ...
    'Color','w','XColor','k','YColor','k','XTick',[],'YTick',[],'Box','on');
s.ax(2) = axes('Units','normalized','Position',[0.347 axY axW axH], ...
    'Color','w','XColor','k','YColor','k','XTick',[],'YTick',[],'Box','on');
s.ax(3) = axes('Units','normalized','Position',[0.684 axY axW axH], ...
    'Color','w','XColor','k','YColor','k','XTick',[],'YTick',[],'Box','on');
for k = 1:3
    title(s.ax(k), s.camLabels{k}, 'Color','k','FontSize',10,'FontWeight','bold');
end

%% ── Flag overview strip ──────────────────────────────────────────────────
s.hStrip = axes('Units','normalized','Position',[0.010 0.100 0.980 0.055], ...
    'Color',[0.97 0.97 0.97],'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4], ...
    'XTick',[],'YTick',[],'Box','on','XLim',[0 1],'YLim',[0 1]);
title(s.hStrip, ...
    'Frame overview  —  flagged = red  |  current = blue  |  click to jump', ...
    'Color',[0.3 0.3 0.3],'FontSize',7.5);
set(s.hStrip,'ButtonDownFcn',@stripClicked);

%% ── Bottom controls ──────────────────────────────────────────────────────
yR1 = 54;   % upper bottom row  (navigation)
yR2 = 12;   % lower bottom row  (flag controls)

% Navigation row
uicontrol('Style','pushbutton','String','◀◀','Position',[8   yR1 42 28], ...
    'BackgroundColor',BTN,'ForegroundColor',FG,'FontSize',9, ...
    'TooltipString','Jump back 10 frames  (Shift+←)', ...
    'Callback',@(~,~) step(-10));
uicontrol('Style','pushbutton','String','◀', 'Position',[54  yR1 42 28], ...
    'BackgroundColor',BTN,'ForegroundColor',FG,'FontSize',9, ...
    'TooltipString','Previous frame  (←)', ...
    'Callback',@(~,~) step(-1));
uicontrol('Style','pushbutton','String','▶', 'Position',[100 yR1 42 28], ...
    'BackgroundColor',BTN,'ForegroundColor',FG,'FontSize',9, ...
    'TooltipString','Next frame  (→)', ...
    'Callback',@(~,~) step(+1));
uicontrol('Style','pushbutton','String','▶▶','Position',[146 yR1 42 28], ...
    'BackgroundColor',BTN,'ForegroundColor',FG,'FontSize',9, ...
    'TooltipString','Jump forward 10 frames  (Shift+→)', ...
    'Callback',@(~,~) step(+10));

s.hFrameLabel = uicontrol('Style','text','String','Frame:  — / —', ...
    'Position',[200 yR1+4 290 21],'FontSize',10,'FontWeight','bold', ...
    'HorizontalAlignment','left','BackgroundColor',BG,'ForegroundColor',FG);

s.hSlider = uicontrol('Style','slider', ...
    'Position',[500 yR1+5 670 20], ...
    'Min',1,'Max',2,'Value',1,'SliderStep',[0.01 0.05], ...
    'BackgroundColor',BTN,'Callback',@sliderMoved);

uicontrol('Style','text','String','Zoom ±px:', ...
    'Position',[1185 yR1+4 68 20],'FontSize',8, ...
    'HorizontalAlignment','right','BackgroundColor',BG,'ForegroundColor',[0.4 0.4 0.4]);
s.hZoom = uicontrol('Style','edit','String','90', ...
    'Position',[1258 yR1 58 28],'FontSize',9,'Callback',@zoomChanged);

% Flag controls row
s.hFlagBtn = uicontrol('Style','pushbutton', ...
    'String','  Flag Frame  [Space]', ...
    'Position',[8 yR2 170 32],'FontSize',10,'FontWeight','bold', ...
    'BackgroundColor',BTN,'ForegroundColor',FG, ...
    'Callback',@(~,~) toggleFlag());

uicontrol('Style','pushbutton','String','▲ Prev Flagged', ...
    'Position',[190 yR2 130 32],'FontSize',9, ...
    'BackgroundColor',BTN,'ForegroundColor',FG, ...
    'TooltipString','Jump to previous flagged frame  (↑)', ...
    'Callback',@(~,~) jumpFlagged(-1));
uicontrol('Style','pushbutton','String','Next Flagged ▼', ...
    'Position',[330 yR2 130 32],'FontSize',9, ...
    'BackgroundColor',BTN,'ForegroundColor',FG, ...
    'TooltipString','Jump to next flagged frame  (↓)', ...
    'Callback',@(~,~) jumpFlagged(+1));

s.hFlagCount = uicontrol('Style','text','String','Flagged:  0 / 0', ...
    'Position',[475 yR2+6 200 22],'FontSize',10,'FontWeight','bold', ...
    'HorizontalAlignment','left','BackgroundColor',BG,'ForegroundColor',[0.65 0.35 0.0]);

uicontrol('Style','pushbutton','String','Save Flags  [S]', ...
    'Position',[690 yR2 150 32],'FontSize',9,'FontWeight','bold', ...
    'BackgroundColor',[0.22 0.40 0.22],'ForegroundColor',FG, ...
    'Callback',@saveFlags);
uicontrol('Style','pushbutton','String','Load Flags', ...
    'Position',[850 yR2 120 32],'FontSize',9, ...
    'BackgroundColor',BTN,'ForegroundColor',FG,'Callback',@loadFlags);
uicontrol('Style','pushbutton','String','Clear All', ...
    'Position',[980 yR2 100 32],'FontSize',9, ...
    'BackgroundColor',[0.40 0.22 0.22],'ForegroundColor',FG,'Callback',@clearFlags);

% Vector toggles
s.hVecBtn = uicontrol('Style','togglebutton','String','Vectors [V]', ...
    'Position',[1090 yR2 105 32],'FontSize',9,'Value',1, ...
    'BackgroundColor',BTN,'ForegroundColor',FG,'Callback',@toggleVectors);
s.hChordBtn = uicontrol('Style','togglebutton','String','Chord [C]', ...
    'Position',[1200 yR2 90 32],'FontSize',9,'Value',0, ...
    'BackgroundColor',BTN,'ForegroundColor',FG,'Callback',@toggleChord);

% Colour legend
uicontrol('Style','text','String','— A-hat (yel)  — span (R/B)  — chord (org/cyan)', ...
    'Position',[1300 yR2+6 220 20],'FontSize',7.5, ...
    'HorizontalAlignment','left','BackgroundColor',BG,'ForegroundColor',[0.25 0.25 0.25]);

guidata(fig, s);
if ~isempty(resultsPath), doLoad(resultsPath); end

%% ════════════════════════════════════════════════════════════════════════
%%  CALLBACKS
%% ════════════════════════════════════════════════════════════════════════

    function loadFile(~,~)
        [fn, fp] = uigetfile('*_results.mat', 'Select *_results.mat file');
        if isequal(fn,0), return; end
        doLoad(fullfile(fp,fn));
    end

    % ── Load & initialise ─────────────────────────────────────────────────
    function doLoad(fpath)
        s = guidata(fig);
        set(s.hFileLabel,'String',['Loading … ' fpath]); drawnow;

        % Load camera calibration + masks from results file
        try
            loaded = load(fpath, 'data','dlt_matrix','all_fly_bw', ...
                          'body_only_bw','params','CM_pos','Nimages');
        catch ME
            errordlg(['Could not load file: ' ME.message], 'Load Error');
            return;
        end

        % ── Try to find a cleaned/corrected data file for res & RESIDX ───
        % The corrected file may have wing-swap corrections applied to RESIDX
        % that aren't in _results.mat.
        % Priority: manually_corrected > iterRefine > cleaned > test > results
        [fdir, fbase] = fileparts(fpath);
        folderBase  = strrep(fbase, '_results', '');   % 'Expr_01_mov_003'
        compactBase = regexprep(folderBase, 'Expr_(\d+)_mov_(\d+)', 'Expr$1mov$2');

        candidates = { ...
            fullfile(fdir, [compactBase '_Data_manually_corrected.mat']), ...
            fullfile(fdir, [folderBase  '_iterRefine.mat']), ...
            fullfile(fdir, [folderBase  '_cleaned.mat']), ...
            fullfile(fdir, [folderBase  '_test.mat']) };

        dataObj    = loaded.data;   % fallback: data from results file
        dataSource = '(results file)';
        for ci = 1:numel(candidates)
            if exist(candidates{ci}, 'file')
                try
                    tmp = load(candidates{ci}, 'data');
                    % handle both top-level 'data' and nested data.data
                    if isfield(tmp, 'data')
                        d = tmp.data ;
                        if isfield(d, 'data'), d = d.data; end
                    else
                        d = tmp ;
                    end
                    % sanity-check that essential fields are present
                    if isfield(d,'res') && isfield(d,'RESIDX')
                        dataObj    = d;
                        [~, dataSource] = fileparts(candidates{ci});
                        dataSource = ['(' dataSource ')'];
                    end
                catch
                end
                break;   % stop at the highest-priority file found
            end
        end

        % ── Vector fields (from cleaned/corrected data) ──────────────────
        s.hasVectors     = false;
        s.bodyCM         = [];   s.rightWingCM   = [];   s.leftWingCM    = [];
        s.AHat           = [];
        s.rightSpanHats  = [];   s.leftSpanHats  = [];
        s.rightChordHats = [];   s.leftChordHats = [];
        if all(isfield(dataObj, {'bodyCM','AHat'}))
            s.bodyCM        = dataObj.bodyCM;
            s.AHat          = dataObj.AHat;
            s.hasVectors    = true;
            if isfield(dataObj,'rightWingCM'),    s.rightWingCM    = dataObj.rightWingCM;    end
            if isfield(dataObj,'leftWingCM'),     s.leftWingCM     = dataObj.leftWingCM;     end
            if isfield(dataObj,'rightSpanHats'),  s.rightSpanHats  = dataObj.rightSpanHats;  end
            if isfield(dataObj,'leftSpanHats'),   s.leftSpanHats   = dataObj.leftSpanHats;   end
            if isfield(dataObj,'rightChordHats'), s.rightChordHats = dataObj.rightChordHats; end
            if isfield(dataObj,'leftChordHats'),  s.leftChordHats  = dataObj.leftChordHats;  end
        end

        % ── Time axis (fps + startTrackingTime) ──────────────────────────
        if isfield(loaded.params,'fps'),               s.fps = loaded.params.fps;                       end
        if isfield(loaded.params,'startTrackingTime'), s.startTrackingTime = loaded.params.startTrackingTime; end

        % Default starting frame: t = -20 ms (clamped to valid range)
        if ~isnan(s.fps) && ~isnan(s.startTrackingTime)
            absTarget = round(-0.020 * s.fps);          % absolute frame at -20 ms
            fTarget   = absTarget - s.startTrackingTime + 1;
            s.frameCurr = max(1, min(loaded.Nimages, fTarget));
        end

        s.loaded      = true;
        s.resultsPath = fpath;
        s.data        = dataObj;
        s.dlt_matrix   = loaded.dlt_matrix;
        s.all_fly_bw   = loaded.all_fly_bw;
        s.body_only_bw = loaded.body_only_bw;
        s.params      = loaded.params;
        s.CM_pos      = loaded.CM_pos;
        s.Nimages     = loaded.Nimages;
        s.frameCurr   = 1;
        s.flagged     = false(loaded.Nimages, 1);

        % Image dimensions
        dpx = loaded.params.detectorLengthPix;
        s.imH = dpx(1);
        s.imW = dpx(end);   % handles scalar or [H W]

        % Camera indices from params
        s.camIdx = [loaded.params.YZ, loaded.params.XZ, loaded.params.XY];
        % DLT columns stay [2 1 3] for the default order=[2 1 3] convention

        % Precompute per-frame row ranges into data.res  (O(M) total)
        [s.startIdx, s.endIdx] = buildFrameIndex( ...
            double(dataObj.res(:,1)), ...
            double(loaded.params.firstTrackableFrame), ...
            loaded.Nimages);

        % Slider range
        N = loaded.Nimages;
        set(s.hSlider, 'Min',1,'Max',max(N,2),'Value',1, ...
            'SliderStep',[1/max(N-1,1), min(10/max(N-1,1),1)]);

        if s.hasVectors
            vecStr = sprintf('vectors: YES (%d frames of AHat)', size(s.AHat,1));
        else
            vecStr = 'vectors: NO  (no cleaned/test/iterRefine file found, or missing fields)';
        end
        set(s.hFileLabel,'String', ...
            sprintf('%s   voxels: %s   |   %s', fpath, dataSource, vecStr));
        guidata(fig, s);
        refreshDisplay();
    end

    % ── Frame navigation ──────────────────────────────────────────────────
    function sliderMoved(~,~)
        s = guidata(fig);
        if ~s.loaded, return; end
        s.frameCurr = max(1, min(s.Nimages, round(get(s.hSlider,'Value'))));
        guidata(fig, s);
        refreshDisplay();
    end

    function step(delta)
        s = guidata(fig);
        if ~s.loaded, return; end
        s.frameCurr = max(1, min(s.Nimages, s.frameCurr + delta));
        guidata(fig, s);
        refreshDisplay();
    end

    function jumpFlagged(dir)
        s = guidata(fig);
        if ~s.loaded, return; end
        if dir > 0
            later = find(s.flagged(s.frameCurr+1 : end));
            if ~isempty(later), s.frameCurr = s.frameCurr + later(1); end
        else
            earlier = find(s.flagged(1 : s.frameCurr-1));
            if ~isempty(earlier), s.frameCurr = earlier(end); end
        end
        guidata(fig, s);
        refreshDisplay();
    end

    function stripClicked(~, event)
        % Click on the overview strip to jump to that frame
        s = guidata(fig);
        if ~s.loaded, return; end
        xFrac = event.IntersectionPoint(1);          % fraction 0–1
        s.frameCurr = max(1, min(s.Nimages, round(1 + xFrac * (s.Nimages-1))));
        guidata(fig, s);
        refreshDisplay();
    end

    % ── Flagging ──────────────────────────────────────────────────────────
    function toggleFlag()
        s = guidata(fig);
        if ~s.loaded, return; end
        s.flagged(s.frameCurr) = ~s.flagged(s.frameCurr);
        guidata(fig, s);
        refreshDisplay();
    end

    function clearFlags(~,~)
        s = guidata(fig);
        if ~s.loaded, return; end
        if strcmp(questdlg('Clear all flags?','Confirm','Yes','No','No'),'Yes')
            s.flagged(:) = false;
            guidata(fig, s);
            refreshDisplay();
        end
    end

    % ── Save / Load flags ─────────────────────────────────────────────────
    function saveFlags(~,~)
        s = guidata(fig);
        if ~s.loaded, return; end
        badFramesList = find(s.flagged); %#ok<NASGU>
        [fdir, fbase] = fileparts(s.resultsPath);
        base = strrep(fbase, '_results', '');
        defPath = fullfile(fdir, [base '_badFrames.mat']);
        [fn, fp] = uiputfile('*.mat', 'Save bad frames list', defPath);
        if isequal(fn,0), return; end
        save(fullfile(fp,fn), 'badFramesList');
        fprintf('Saved %d flagged frames → %s\n', numel(badFramesList), fullfile(fp,fn));
        title(s.hStrip, sprintf( ...
            'Saved %d flagged frames → %s', numel(badFramesList), fn), ...
            'Color',[0.0 0.55 0.0],'FontSize',7.5);
    end

    function loadFlags(~,~)
        s = guidata(fig);
        if ~s.loaded, errordlg('Load a results file first.'); return; end
        [fn, fp] = uigetfile('*_badFrames.mat','Load bad frames list');
        if isequal(fn,0), return; end
        d = load(fullfile(fp,fn),'badFramesList');
        s.flagged(:) = false;
        idx = d.badFramesList;
        s.flagged(idx(idx>=1 & idx<=s.Nimages)) = true;
        guidata(fig, s);
        refreshDisplay();
    end

    % ── Zoom edit ─────────────────────────────────────────────────────────
    function zoomChanged(src,~)
        s = guidata(fig);
        v = str2double(get(src,'String'));
        if ~isnan(v) && v > 0
            s.zoom = round(v);
            guidata(fig, s);
            refreshDisplay();
        end
    end

    % ── Keyboard shortcuts ────────────────────────────────────────────────
    function keyPressed(~, event)
        switch event.Key
            case 'rightarrow'
                if any(strcmp(event.Modifier,'shift'))
                    step(+10);
                else
                    step(+1);
                end
            case 'leftarrow'
                if any(strcmp(event.Modifier,'shift'))
                    step(-10);
                else
                    step(-1);
                end
            case 'uparrow',   jumpFlagged(-1);
            case 'downarrow', jumpFlagged(+1);
            case 'space',     toggleFlag();
            case 's',         saveFlags([],[]);
            case 'v'
                s = guidata(fig);
                s.showVectors = ~s.showVectors;
                set(s.hVecBtn,'Value', s.showVectors);
                guidata(fig, s);
                refreshDisplay();
            case 'c'
                s = guidata(fig);
                s.showChord = ~s.showChord;
                set(s.hChordBtn,'Value', s.showChord);
                guidata(fig, s);
                refreshDisplay();
        end
    end

    function toggleVectors(src,~)
        s = guidata(fig);
        s.showVectors = logical(get(src,'Value'));
        guidata(fig, s);
        refreshDisplay();
    end

    function toggleChord(src,~)
        s = guidata(fig);
        s.showChord = logical(get(src,'Value'));
        guidata(fig, s);
        refreshDisplay();
    end

%% ════════════════════════════════════════════════════════════════════════
%%  DISPLAY
%% ════════════════════════════════════════════════════════════════════════

    function refreshDisplay()
        s = guidata(fig);
        if ~s.loaded, return; end
        f = s.frameCurr;
        N = s.Nimages;

        % ── Update text / slider ─────────────────────────────────────────
        if ~isnan(s.fps) && ~isnan(s.startTrackingTime)
            absFrame = f - 1 + s.startTrackingTime;
            t_ms     = absFrame / s.fps * 1000;
            set(s.hFrameLabel, 'String', ...
                sprintf('Frame:  %d / %d   t = %+6.2f ms', f, N, t_ms));
        else
            set(s.hFrameLabel, 'String', sprintf('Frame:  %d / %d', f, N));
        end
        set(s.hSlider, 'Value', f);
        set(s.hFlagCount, 'String', sprintf('Flagged:  %d / %d', sum(s.flagged), N));

        if s.flagged(f)
            set(s.hFlagBtn, 'String','  ✓ Unflag  [Space]', ...
                'BackgroundColor',[0.60 0.18 0.18]);
        else
            set(s.hFlagBtn, 'String','  Flag Frame  [Space]', ...
                'BackgroundColor',[0.85 0.85 0.85],'ForegroundColor',[0.10 0.10 0.10]);
        end

        % ── Get voxels for this frame ─────────────────────────────────────
        i1 = s.startIdx(f);
        i2 = s.endIdx(f);
        if i1 > 0 && i2 >= i1
            vox   = double(s.data.res(i1:i2, 2:4));
            world = vox * s.params.voxelSize;           % Nx3 in metres
            ridx  = s.data.RESIDX(i1:i2, :);
            iB = logical(ridx(:, s.data.bodyInd));
            iR = logical(ridx(:, s.data.rightWingInd));
            iL = logical(ridx(:, s.data.leftWingInd));
        else
            world = zeros(0,3);
            iB = false(0,1); iR = iB; iL = iB;
        end

        % ── Draw each camera view ─────────────────────────────────────────
        imH = s.imH;  imW = s.imW;

        for ic = 1:3
            axh    = s.ax(ic);
            camI   = s.camIdx(ic);
            dltCol = s.dltCols(ic);

            % Composite background:
            %   body region  (green) = body_only_bw
            %   wing region  (grey)  = all_fly_bw minus body
            %   background   (black) = outside fly
            flyMask  = getImage4D(s.all_fly_bw,   camI, f);
            bodyMask = getImage4D(s.body_only_bw,  camI, f);
            wingMask = flyMask & ~bodyMask;
            R = uint8(bodyMask) * 25  + uint8(wingMask) * 55;
            G = uint8(bodyMask) * 100 + uint8(wingMask) * 55;
            B = uint8(bodyMask) * 25  + uint8(wingMask) * 55;

            % Project and paint voxels
            if ~isempty(world)
                pix = dlt_inverse(s.dlt_matrix(:, dltCol), world);
                pu  = round(pix(:,1));                 % pixel column
                pv  = round(imH - pix(:,2) + 1);      % pixel row (flip y)
                ok  = pu>=1 & pu<=imW & pv>=1 & pv<=imH;

                % Alpha-blend voxels over background so body/wing mask stays visible
                % body=white, right=red, left=blue
                alpha = 0.55;
                if any(ok & iB)
                    idx = sub2ind([imH imW], pv(ok&iB), pu(ok&iB));
                    R(idx) = uint8(double(R(idx))*(1-alpha) + 210*alpha);
                    G(idx) = uint8(double(G(idx))*(1-alpha) + 210*alpha);
                    B(idx) = uint8(double(B(idx))*(1-alpha) + 210*alpha);
                end
                if any(ok & iR)
                    idx = sub2ind([imH imW], pv(ok&iR), pu(ok&iR));
                    R(idx) = uint8(double(R(idx))*(1-alpha) + 255*alpha);
                    G(idx) = uint8(double(G(idx))*(1-alpha) + 55*alpha);
                    B(idx) = uint8(double(B(idx))*(1-alpha) + 55*alpha);
                end
                if any(ok & iL)
                    idx = sub2ind([imH imW], pv(ok&iL), pu(ok&iL));
                    R(idx) = uint8(double(R(idx))*(1-alpha) + 55*alpha);
                    G(idx) = uint8(double(G(idx))*(1-alpha) + 100*alpha);
                    B(idx) = uint8(double(B(idx))*(1-alpha) + 255*alpha);
                end
            end

            imshow(cat(3, R, G, B), 'Parent', axh);

            % ── Vector overlay (A-hat, span, chord) ─────────────────────
            if s.showVectors && s.hasVectors
                hold(axh,'on');
                % body long axis
                drawVec(axh, s.bodyCM, s.AHat, s.vecScaleBody, ...
                        [1.0 0.95 0.30], 2.0, dltCol, imH, imW, f, s);
                % right wing span
                if ~isempty(s.rightSpanHats) && ~isempty(s.rightWingCM)
                    drawVec(axh, s.rightWingCM, s.rightSpanHats, s.vecScaleSpan, ...
                            [1.0 0.30 0.30], 1.7, dltCol, imH, imW, f, s);
                end
                % left wing span
                if ~isempty(s.leftSpanHats) && ~isempty(s.leftWingCM)
                    drawVec(axh, s.leftWingCM, s.leftSpanHats, s.vecScaleSpan, ...
                            [0.30 0.60 1.0], 1.7, dltCol, imH, imW, f, s);
                end
                % chord vectors (only when chord toggle is on)
                if s.showChord
                    if ~isempty(s.rightChordHats) && ~isempty(s.rightWingCM)
                        drawVec(axh, s.rightWingCM, s.rightChordHats, s.vecScaleChord, ...
                                [1.0 0.55 0.10], 1.4, dltCol, imH, imW, f, s);
                    end
                    if ~isempty(s.leftChordHats) && ~isempty(s.leftWingCM)
                        drawVec(axh, s.leftWingCM, s.leftChordHats, s.vecScaleChord, ...
                                [0.20 0.85 0.85], 1.4, dltCol, imH, imW, f, s);
                    end
                end
                hold(axh,'off');
            end

            % Zoom to fly using CM_pos
            if ~isempty(s.CM_pos) && ~isnan(s.CM_pos(camI, f, 1))
                xc = s.CM_pos(camI, f, 1);
                yc = s.CM_pos(camI, f, 2);
                hw = s.zoom;
                set(axh, 'XLim',[xc-hw, xc+hw], 'YLim',[yc-hw, yc+hw]);
            end

            % Red border if flagged, white otherwise
            if s.flagged(f)
                set(axh,'XColor',[0.9 0.2 0.2],'YColor',[0.9 0.2 0.2],'LineWidth',2.5);
            else
                set(axh,'XColor','k','YColor','k','LineWidth',0.5);
            end
            title(axh, s.camLabels{ic}, 'Color','k','FontSize',10,'FontWeight','bold');
        end

        updateStrip(s);
    end

    % ── Flag overview strip ───────────────────────────────────────────────
    function updateStrip(s)
        axh = s.hStrip;
        cla(axh);
        N = s.Nimages;
        if N < 2, return; end

        set(axh, 'XLim',[0 1],'YLim',[0 1], ...
            'Color',[0.97 0.97 0.97],'XTick',[],'YTick',[], ...
            'ButtonDownFcn',@stripClicked);
        hold(axh,'on');

        % Draw flagged frames as red patches
        flagInds = find(s.flagged);
        if ~isempty(flagInds)
            x0 = (flagInds - 1) / (N - 1);
            dx = 1 / (N - 1);
            for k = 1:numel(flagInds)
                patch(axh, x0(k)+[0 dx dx 0], [0 0 1 1], ...
                    [0.85 0.22 0.22], 'EdgeColor','none');
            end
        end

        % Current frame as cyan line
        xc = (s.frameCurr - 1) / (N - 1);
        plot(axh, [xc xc], [0 1], '-', 'Color',[0.05 0.40 0.85], 'LineWidth', 2);

        hold(axh,'off');
        title(axh, ...
            sprintf('Frame overview  —  flagged = %d red  |  current = cyan  |  click to jump', sum(s.flagged)), ...
            'Color',[0.3 0.3 0.3],'FontSize',7.5);
    end

%% ════════════════════════════════════════════════════════════════════════
%%  HELPERS
%% ════════════════════════════════════════════════════════════════════════

    function drawVec(axh, cmMat, hatMat, scale, color, lw, dltCol, imH, imW, f, s)
        % Project a single 3-D vector (cm → cm+hat*scale) into camera view and draw it.
        % cmMat is in VOXEL units (matching data.res), hatMat is unit vector, scale is metres.
        if size(cmMat,1) < f || size(hatMat,1) < f, return; end
        p0 = cmMat(f, :) * s.params.voxelSize;   % voxel → metres
        v  = hatMat(f, :);
        if any(~isfinite(p0)) || any(~isfinite(v)) || all(v == 0), return; end
        % Project both endpoints centred on the CM so the vector reads as a head/tail line.
        pA = p0 - 0.5 * scale * v;
        pB = p0 + 0.5 * scale * v;
        px = dlt_inverse(s.dlt_matrix(:, dltCol), [pA; pB]);
        % Convert to image-row coords (flip y)
        u  = px(:,1);
        v2 = imH - px(:,2) + 1;
        % Clip to image bounds for safety (still allow lines that exit zoom window)
        if all(u >= -imW) && all(u <= 2*imW) && all(v2 >= -imH) && all(v2 <= 2*imH)
            plot(axh, u, v2, '-', 'Color', color, 'LineWidth', lw);
            plot(axh, u(2), v2(2), '.', 'Color', color, 'MarkerSize', 9);
        end
    end

    function [si, ei] = buildFrameIndex(absVec, firstAbs, N)
        % Map each 1-based frame index → [start_row, end_row] in data.res.
        % Runs in O(M) via unique().
        si = zeros(N, 1, 'int32');
        ei = zeros(N, 1, 'int32');
        [uFrames, ia] = unique(absVec, 'first');
        [~,        iz] = unique(absVec, 'last');
        for j = 1:numel(uFrames)
            k = uFrames(j) - firstAbs + 1;
            if k >= 1 && k <= N
                si(k) = ia(j);
                ei(k) = iz(j);
            end
        end
    end

end % hullInspectionGUI
