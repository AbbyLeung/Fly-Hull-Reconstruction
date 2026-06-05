classdef correctionGUI_sam_App_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        figure1                 matlab.ui.Figure
        uitoolbar1              matlab.ui.container.Toolbar
        ui_zoom_in              matlab.ui.container.toolbar.ToggleTool
        ui_zoom_out             matlab.ui.container.toolbar.ToggleTool
        ui_rotate               matlab.ui.container.toolbar.ToggleTool
        ui_pan                  matlab.ui.container.toolbar.ToggleTool
        ui_save                 matlab.ui.container.toolbar.PushTool
        roll_view               matlab.ui.control.Button
        clear_data              matlab.ui.control.Button
        ffwd                    matlab.ui.control.Button
        fwd                     matlab.ui.control.Button
        back                    matlab.ui.control.Button
        bback                   matlab.ui.control.Button
        load_data               matlab.ui.control.Button
        open_data_dir           matlab.ui.control.Button
        data_dir                matlab.ui.control.ListBox
        correction_panel        matlab.ui.container.Panel
        r_eta_view              matlab.ui.control.Button
        phi_view                matlab.ui.control.Button
        l_theta_view            matlab.ui.control.Button
        phi_view_2              matlab.ui.control.Button
        l_eta_view              matlab.ui.control.Button
        r_theta_view            matlab.ui.control.Button
        EZViewButton            matlab.ui.control.Button
        l_eta_inc               matlab.ui.control.Button
        l_theta_inc             matlab.ui.control.Button
        l_phi_inc               matlab.ui.control.Button
        l_z_inc                 matlab.ui.control.Button
        l_y_inc                 matlab.ui.control.Button
        l_x_inc                 matlab.ui.control.Button
        r_eta_inc               matlab.ui.control.Button
        r_theta_inc             matlab.ui.control.Button
        r_phi_inc               matlab.ui.control.Button
        r_z_inc                 matlab.ui.control.Button
        r_y_inc                 matlab.ui.control.Button
        rx_inc                  matlab.ui.control.Button
        body_roll_inc           matlab.ui.control.Button
        body_pitch_inc          matlab.ui.control.Button
        body_phi_inc            matlab.ui.control.Button
        body_z_inc              matlab.ui.control.Button
        l_eta_dec               matlab.ui.control.Button
        l_theta_dec             matlab.ui.control.Button
        l_phi_dec               matlab.ui.control.Button
        l_z_dec                 matlab.ui.control.Button
        l_y_dec                 matlab.ui.control.Button
        l_x_dec                 matlab.ui.control.Button
        r_eta_dec               matlab.ui.control.Button
        r_theta_dec             matlab.ui.control.Button
        r_phi_dec               matlab.ui.control.Button
        r_z_dec                 matlab.ui.control.Button
        r_y_dec                 matlab.ui.control.Button
        rx_dec                  matlab.ui.control.Button
        body_roll_dec           matlab.ui.control.Button
        body_pitch_dec          matlab.ui.control.Button
        body_phi_dec            matlab.ui.control.Button
        body_z_dec              matlab.ui.control.Button
        body_y_inc              matlab.ui.control.Button
        body_y_dec              matlab.ui.control.Button
        body_x_inc              matlab.ui.control.Button
        body_x_dec              matlab.ui.control.Button
        save_roll               matlab.ui.control.StateButton
        save_manualCorrRangeMS  matlab.ui.control.Button
        cluster_wings           matlab.ui.control.Button
        l_z_label               matlab.ui.control.Label
        l_y_label               matlab.ui.control.Label
        l_x_label               matlab.ui.control.Label
        r_z_label               matlab.ui.control.Label
        r_y_label               matlab.ui.control.Label
        r_x_label               matlab.ui.control.Label
        roll_label              matlab.ui.control.Label
        pitch_label             matlab.ui.control.Label
        yaw_label               matlab.ui.control.Label
        z_cm_label              matlab.ui.control.Label
        y_cm_label              matlab.ui.control.Label
        swap_wings              matlab.ui.control.Button
        ignoreFrame             matlab.ui.control.StateButton
        body_roll               matlab.ui.control.Slider
        body_pitch              matlab.ui.control.Slider
        body_phi                matlab.ui.control.Slider
        left_wing_label         matlab.ui.control.Label
        right_wing_label        matlab.ui.control.Label
        l_eta                   matlab.ui.control.Slider
        l_theta                 matlab.ui.control.Slider
        l_phi                   matlab.ui.control.Slider
        l_z                     matlab.ui.control.Slider
        l_y                     matlab.ui.control.Slider
        l_x                     matlab.ui.control.Slider
        r_eta                   matlab.ui.control.Slider
        r_theta                 matlab.ui.control.Slider
        r_phi                   matlab.ui.control.Slider
        r_z                     matlab.ui.control.Slider
        r_y                     matlab.ui.control.Slider
        rx                      matlab.ui.control.Slider
        body_z                  matlab.ui.control.Slider
        body_y                  matlab.ui.control.Slider
        x_cm_label              matlab.ui.control.Label
        body_label              matlab.ui.control.Label
        body_x                  matlab.ui.control.Slider
        frame_info              matlab.ui.control.Label
        main_axes               matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        EZViewApp % Description
        frameCurr
        tmsCurr
    end
    
    methods (Access = private)
        function clearDisplay(app, hObject)
            % -------------------------------------
            % remove current fly plot from handles
            
            handles = guidata(hObject);
            
            x_data = [NaN, NaN] ;
            y_data = [NaN, NaN] ;
            z_data = [NaN, NaN] ;
            
            % voxels
            set(handles.hBody, 'XData', NaN, 'YData', NaN,'ZData', NaN);
            set(handles.hRightWing, 'XData', NaN,'YData', NaN,'ZData', NaN);
            set(handles.hLeftWing, 'XData', NaN,'YData', NaN, 'ZData', NaN);
            
            % body vectors
            set(handles.hAhat, 'XData', x_data,'YData', y_data,'ZData', z_data);
            set(handles.hStubline, 'XData', x_data,'YData', y_data,'ZData', z_data);
            set(handles.hNr, 'XData', x_data,'YData', y_data,'ZData', z_data);
            set(handles.hPh, 'XData', x_data,'YData', y_data,'ZData', z_data);
            
            % wing vectors
            set(handles.hSpanRight, 'XData', x_data,'YData', y_data,'ZData', z_data);
            set(handles.hChordRight, 'XData', x_data,'YData', y_data,'ZData', z_data);
            set(handles.hSpanLeft, 'XData', x_data,'YData', y_data,'ZData', z_data);
            set(handles.hChordLeft, 'XData', x_data,'YData', y_data, 'ZData', z_data);
            
            set(handles.frame_info,'String','Frame Info')
            
            guidata(hObject, handles);
        end
        
        function enableButtons(app, hObject)
            % ----------------------------------------------
            % turn buttons on
            
            handles = guidata(hObject) ;
            %slider_handles = findall(0,'Style','Slider') ;
            slider_handles = [handles.body_roll, handles.body_pitch, handles.body_phi, handles.l_eta, handles.l_theta, handles.l_phi, handles.l_z,...
                handles.l_y, handles.l_x, handles.r_eta, handles.r_theta, handles.r_phi, handles.r_z, handles.rx, handles.body_z,...
                handles.body_y, handles.body_x, handles.r_y];
            set(slider_handles,'Enable','on')
            %button_handles = findall(0,'Style','pushbutton') ;
            button_handles = [handles.clear_data, handles.ffwd, handles.fwd, handles.back, handles.bback, handles.load_data, handles.open_data_dir,...
                handles.save_manualCorrRangeMS, handles.cluster_wings, handles.swap_wings, app.body_x_dec, app.body_x_inc, app.body_y_dec, app.body_y_inc,...
                app.body_z_dec, app.body_phi_dec, app.body_pitch_dec, app.body_roll_dec, app.rx_dec, app.r_y_dec, app.r_z_dec, app.r_phi_dec,...
                app.r_theta_dec, app.r_eta_dec, app.l_x_dec, app.l_y_dec, app.l_z_dec, app.l_phi_dec, app.l_theta_dec, app.l_eta_dec, app.body_z_inc,...
                app.body_phi_inc, app.body_pitch_inc, app.body_roll_inc, app.rx_inc, app.r_y_inc, app.r_z_inc, app.r_phi_inc, app.r_theta_inc, app.r_eta_inc,...
                app.l_x_inc, app.l_y_inc, app.l_z_inc, app.l_phi_inc, app.l_theta_inc, app.l_eta_inc, app.EZViewButton, app.phi_view, app.r_theta_view, ...
                app.l_theta_view, app.r_eta_view, app.l_eta_view, app.roll_view];
            set(button_handles,'Enable','on')
            %toggle_button_handles = findall(0,'Style','togglebutton') ;
            %set(toggle_button_handles,'Enable','on')
            state_button_handles = [handles.save_roll,handles.ignoreFrame];
            set(state_button_handles,'Enable','on')
        end
        
        function handles = initPlots(app, handles)
            % ----------------------------------
            % initialize plot handles
            
            % initialize handle fields for things to plot
            
            scale = 4 ; % 2; % 4
            handles.scale = scale ;
            handles.Lstub = 3.0*scale ;
            handles.Lbar  = 10*scale ;
            
            %handles.ignoreFrames = [] ;
            
            % flags to mark whether or not things need to be updated
            handles.bodyChangeFlag = false ;
            handles.rightWingChangeFlag = false ;
            handles.leftWingChangeFlag = false ;
            
            handles.AhatChangeFlag = false ;
            handles.stublineChangeFlag = false ;
            handles.nrChangeFlag = false ;
            handles.phChangeFlag = false ;
            
            handles.spanRightChangeFlag = false ;
            handles.chordRightChangeFlag = false ;
            handles.spanLeftChangeFlag = false ;
            handles.chordLeftChangeFlag = false ;
            
            % intialize plot elements
            
            hold(handles.main_axes, 'on')
            handles.hBody = plot3(handles.main_axes, NaN, NaN, NaN, 'g.',...
                'MarkerSize',handles.marker_size) ;
            handles.hRightWing = plot3(handles.main_axes, NaN, NaN, NaN, 'r.',...
                'MarkerSize',handles.marker_size) ;
            handles.hLeftWing = plot3(handles.main_axes, NaN, NaN, NaN, 'b.',...
                'MarkerSize',handles.marker_size) ;
            
            handles.hAhat = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
                'Color','r','LineWidth',8);
            handles.hStubline = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
                'Color','k','LineWidth',8);
            handles.hNr = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
                'Color','k','LineWidth',8);
            handles.hPh = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
                'Color','r','LineWidth',8);
            
            handles.hSpanRight = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
                'Color','k','LineWidth',4);
            handles.hChordRight = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
                'Color','b','LineWidth',4);
            handles.hSpanLeft = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
                'Color','k','LineWidth',4);
            handles.hChordLeft = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
                'Color','r','LineWidth',4);
            
            hold(handles.main_axes, 'off')
            box(handles.main_axes, 'on')
            grid(handles.main_axes, 'on')
            axis(handles.main_axes, 'equal')
            rotate3d on
            
            %initialize view
            azview = -43 ;
            elview =  26 ;
            view(azview,elview) ;
        end
        
        function handles = loadData(app, hObject, handles)
            % ---------------------------------------------------------------------
            % function to load in a fly analysis data file
            
            if ~isempty(handles.datapath_curr)
            
                % load data file
                data = importdata(handles.datapath_curr) ;
                if isfield(data,'data')
                    data = data.data ;
                end
            
                % do we switch to body frame coordinates?
                if handles.bodyFrameFlag
                   data = labToBodyFrame(data, handles.largePertFlag) ;
                   handles.bodyFrameRotMats = data.bodyFrameRotMats ;
                   handles.bodyCM_orig = data.bodyCM_orig ;
                end
                %================================================
                % add fields to handles for relevant plotting data
                %================================================
                % vectors for body/wing orientation
                handles.chord1Hats = data.rightChordHats ;
                handles.chord2Hats = data.leftChordHats ;
                handles.span1Hats  = data.rightSpanHats ;
                handles.span2Hats  = data.leftSpanHats ;
                handles.rollHats   = data.AHat ; % long body axis
            
                handles.Nimages     = data.Nimages ;
                handles.wingLength  = data.wingLength ;
                psiVec  = atan2(data.AHat(:,2),data.AHat(:,1)); % body angle with respect to x axis
                handles.psiHats  = [-sin(psiVec)  cos(psiVec)  zeros(data.Nimages,1)];
            
                phi1Vec   = atan2(data.rightSpanHats(:,2),data.rightSpanHats(:,1));
                handles.phi1Hats = [-sin(phi1Vec) cos(phi1Vec) zeros(data.Nimages,1)];
            
                phi2Vec   = atan2(data.leftSpanHats(:,2),data.leftSpanHats(:,1));
                handles.phi2Hats = [-sin(phi2Vec) cos(phi2Vec) zeros(data.Nimages,1)];
            
                % time information
                df = diff(data.res(:,1)) ;
                frameStartInd = [1 ; find(df==1)+1] ;
                handles.frameStartInd = frameStartInd ;
                handles.frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;
                handles.noffset       = data.params.startTrackingTime - data.params.firstTrackableFrame  ;
                handles.params        =data.params ;
            
                tvec = (0:data.Nimages-1) + data.params.startTrackingTime ;
                handles.tvec = tvec / 8000 * 1000 ; % in ms
            
                app.tmsCurr = handles.tmsCurr_default ;
            
                % voxel info
                handles.res = data.res ;
                handles.RESIDX = data.RESIDX ;
            
                % cm info
                handles.bodyCM = data.bodyCM ;
                handles.rightWingCM = data.rightWingCM ;
                handles.leftWingCM = data.leftWingCM ;
            
                % wing tip info
                handles.rightWingTips = data.rightWingTips ;
                handles.leftWingTips = data.leftWingTips ;
            
                % indexing for body/wings
                handles.bodyInd = data.bodyInd ;
                handles.rightWingInd = data.rightWingInd ;
                handles.leftWingInd = data.leftWingInd ;
            
                %================================================
                % mark that all plot objects need to be updated
                %================================================
                handles = updateChangeFlags(app, handles,'all') ;
            
                % mark that data is currently loaded
                handles.dataFlag = true ;
                datapath_curr = handles.datapath_curr ;
                [filepath, name, ext] = fileparts(datapath_curr) ;
                if contains(name,'manually_corrected')
                    handles.output_path = datapath_curr ;
                else
                    name_split = strsplit(name, '_') ;
                    ExprNumStr = name_split{2} ;
                    MovNumStr = name_split{4} ;
                    name_new = ['Expr' ExprNumStr 'mov' MovNumStr ...
                        '_Data_manually_corrected' ext] ;
                    handles.output_path = fullfile(filepath, name_new) ;
                end
            
                % initialize matrices for adjustments/swaps to data
                handles.adjust = zeros(data.Nimages,19);
                handles.rhoFlag   = false(data.Nimages,1);
                handles.ignoreFlag = false(data.Nimages,1) ;
                handles.pitchflipFlag = false(data.Nimages,1) ;
                handles.wingSwapFlag = false(data.Nimages,1) ;
            
                % add in fields from data if they exist
                if isfield(data,'rhoTimes')
                    handles.rhoTimes = data.rhoTimes ;
                    handles.rhoFlag(handles.rhoTimes) = true ;
                else
                    handles.rhoTimes = [] ;
                end
                if isfield(data,'rollVectors')
                    handles.normRolls   = data.rollVectors ;
                else
                    handles.normRolls =  zeros(data.Nimages,3) ;
                end
                if isfield(data,'ignoreFrames')
                   handles.ignoreFrames = data.ignoreFrames ;
                   handles.ignoreFlag(handles.ignoreFrames) = true ;
                else
                    handles.ignoreFrames = [] ;
                end
                % enable sliders so that angles can be adjusted
                enableButtons(app, hObject)
            
                % store data to gui and update display
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            else
                % if there's nothing to load, just print that and do nothing
                disp('No data file selected')
            end
        end
        
        function resetSliders(app, hObject, newFrame)
            % -----------------------------------
            % return sliders to default position
            
            % get handles struct from gui object
            handles = guidata(hObject) ;
            
            % set all sliders back to zero (or where they were from previous correction)
            %slider_handles = findall(0,'Style','Slider') ;
            slider_handles = [handles.body_roll, handles.body_pitch, handles.body_phi, handles.l_eta, handles.l_theta, handles.l_phi, handles.l_z,...
                handles.l_y, handles.l_x, handles.r_eta, handles.r_theta, handles.r_phi, handles.r_z, handles.rx, handles.body_z,...
                handles.body_y, handles.body_x, handles.r_y];
            % set(slider_handles,'Value',0) ;
            for ind = 1:length(slider_handles)
               slider_curr = slider_handles(ind) ;
               tag_curr = get(slider_curr,'Tag') ;
               reset_val = handles.adjust(newFrame, ...
                   handles.adjustTagStruct.(tag_curr)) ;
               set(slider_curr,'Value', reset_val)
               % disp(reset_val)
            end
        end
        
        function resetToggleButtons(app, hObject, newFrame)
            % ----------------------------------------
            % return toggle buttons to default status
            
            % read in handle data
            handles = guidata(hObject) ;
            % get current frame
            frameCurr = newFrame ;
            % find toggle buttons
            %toggle_button_handles = findall(0,'Style','togglebutton') ;
            state_button_handles = [handles.save_roll,handles.ignoreFrame];
            % get default button background color
            defaultColor = get(0,'defaultUicontrolBackgroundColor') ;
            % loop through toggle buttons to update each depending on their value in
            % new frame
            N_buttons = length(state_button_handles) ;
            for ind = 1:N_buttons
                button_handle_curr = state_button_handles(ind) ;
                button_string_curr = get(button_handle_curr,'String') ;
                switch button_string_curr
                    case 'Save Roll'
                        toggleVal = handles.rhoFlag(frameCurr) ;
                        colorStr = 'green' ;
                    case 'Ignore Frame'
                        toggleVal = handles.ignoreFlag(frameCurr) ;
                        colorStr = 'red' ;
                    otherwise
                        toggleVal = 0 ;
                        colorStr = defaultColor ;
                end
            
                % set toggle button value and apply color change
                set(button_handle_curr,'Value',toggleVal) ;
                if toggleVal
                    set(button_handle_curr,'BackgroundColor',colorStr) ;
                else
                    set(button_handle_curr,'BackgroundColor',defaultColor) ;
                end
            end
            guidata(hObject, handles);
        end
        
        function saveData(app, hObject)
            % ---------------------------------
            % save *_manually_corrected file
            
            % function to save manual correction results
            handles = guidata(hObject) ;
            f = waitbar(0,'') ;
            f.Children(end).Title.Interpreter = 'none' ;
            [~, fn_str, ~] = fileparts(handles.output_path) ;
            waitbar(0,f, sprintf('Saving data to %s', fn_str)) ;
            
            % load data file
            data = importdata(handles.datapath_curr) ;
            if isfield(data,'data')
                data = data.data ;
            end
            
            waitbar(0.33)
            % replace values in data structure with adjusted values
            data.bodyCM = handles.bodyCM ;
            data.rightWingCM = handles.rightWingCM ;
            data.leftWingCM  = handles.leftWingCM ;
            
            % vectors
            data.rightChordHats = handles.chord1Hats ;
            data.leftChordHats  = handles.chord2Hats ;
            data.rightSpanHats  = handles.span1Hats ;
            data.leftSpanHats   = handles.span2Hats ;
            data.AHat           = handles.rollHats ;
            
            data.rhoTimes       = handles.rhoTimes ;
            data.rollVectors    = handles.normRolls ;
            data.ignoreFrames   = handles.ignoreFrames ;
            
            data.RESIDX         = handles.RESIDX ;
            
            if isfield(data,'rightWingTips')
                temp1 = data.rightWingTips ;
                data.rightWingTips([handles.wingSwapFlag],:) = ...
                    data.leftWingTips([handles.wingSwapFlag],:) ;
                data.leftWingTips([handles.wingSwapFlag],:) = ...
                    temp1([handles.wingSwapFlag],:) ;
            end
            
            if isfield(handles, 'manualCorrRangeMS')
                data.manualCorrRangeMS = handles.manualCorrRangeMS ;
            end
            
            if handles.bodyFrameFlag
               data.bodyCM_orig = handles.bodyCM_orig ;
               data.bodyFrameRotMats = handles.bodyFrameRotMats ;
               data = bodyToLabFrame(data) ;
            end
            
            % update and close waitbar
            waitbar(0.66)
            save(handles.output_path,'data')
            waitbar(1)
            close(f)
            
            % also calculate angles?
            if handles.calcNewAnglesFlag
                % make new waitbar to alert user of angle calculation
                f_ang = waitbar(0,'') ;
                f_ang.Children(end).Title.Interpreter = 'none' ;
                waitbar(0,f_ang, sprintf('Calculating angles for %s', fn_str)) ;
            
                % calculate angles
                angleSaveFlag = true ;
                anglePlotFlag = true ;
                calcAnglesMain(handles.output_path, handles.largePertFlag, ...
                    angleSaveFlag, anglePlotFlag) ;
            
                % close waitbar
                waitbar(1)
                close(f_ang)
            
                % close plots
                openFigures = findobj('Type','Figure','-not','Tag',...
                    get(handles.output,'Tag'));
                close(openFigures)
            
            end
        end
        
        function [handles, adjust_val] = sliderAdjust(app, hObject, handles)
            % generic function for data adjustment sliders
            tag_curr = get(hObject,'Tag') ;
            slider_val_curr = get(hObject,'Value') ;
            slider_val_prev = ...
                handles.adjust(app.frameCurr,handles.adjustTagStruct.(tag_curr)) ;
            adjust_val = slider_val_curr - slider_val_prev ;
            handles.adjust(app.frameCurr,handles.adjustTagStruct.(tag_curr)) = ...
                slider_val_curr ;
        end
        
        function handles = updateChangeFlags(app, handles, updateType)
            % ---------------------------------------------------------------------
            % function to store whether or not a given fly objet was changed
            
            % flag all variables for update display (useful for frame shifts)
            switch updateType
                case 'all'
                    handles.bodyChangeFlag = true ;
                    handles.rightWingChangeFlag = true ;
                    handles.leftWingChangeFlag = true ;
            
                    handles.AhatChangeFlag = true ;
                    handles.stublineChangeFlag = true ;
                    handles.nrChangeFlag = true ;
                    handles.phChangeFlag = true ;
            
                    handles.spanRightChangeFlag = true ;
                    handles.chordRightChangeFlag = true ;
                    handles.spanLeftChangeFlag = true ;
                    handles.chordLeftChangeFlag = true ;
                case 'bodyvec'
                    handles.AhatChangeFlag = true ;
                    handles.stublineChangeFlag = true ;
                    handles.nrChangeFlag = true ;
                    handles.phChangeFlag = true ;
                case 'rvec'
                    handles.spanRightChangeFlag = true ;
                    handles.chordRightChangeFlag = true ;
                case 'lvec'
                    handles.spanLeftChangeFlag = true ;
                    handles.chordLeftChangeFlag = true ;
                case 'rchord'
                    handles.chordRightChangeFlag = true ;
                case 'lchord'
                    handles.chordLeftChangeFlag = true ;
            end
        end
        
        function updateDisplay(app, hObject)
            % --------------------------------
            % update plot to reflect changes
            
            handles = guidata(hObject);
            
            %============================
            % determine current frame
            %============================
            
                frameCurr = app.frameCurr; 
            if isempty(frameCurr)
                tmsCurr = app.tmsCurr ;
                tvec = handles.tvec ;
                [~, frameCurr] = min(abs(tvec - tmsCurr)) ;
                app.frameCurr = frameCurr ;
                app.tmsCurr = tvec(frameCurr) ;
            end
            
            %
            %============================
            % get voxel info
            %============================
            row1 = handles.frameStartInd(frameCurr) ;
            row2 = handles.frameEndInd(frameCurr) ;
            coords = handles.res(row1:row2,2:4) ; % xyz
            IDX = handles.RESIDX(row1:row2,:) ;  % (isbody, isrightwing, isleftwing)
            bodyRows      = (IDX(:,handles.bodyInd)==1) ;
            rightWingRows = (IDX(:,handles.rightWingInd)==1) ;
            leftWingRows  = (IDX(:,handles.leftWingInd)==1) ;
            
            scale = handles.scale ; % constant
            Lbar = handles.Lbar ; % constant
            
            %============================
            % get body vector info
            %============================
            cb = handles.bodyCM(frameCurr,:) ; % body cm
            rollHat = handles.rollHats(frameCurr,:) ; % body axis unit vector
            
            checksum = sum(handles.normRolls(frameCurr,:));
            if ((~isfinite(checksum) || norm(handles.normRolls(frameCurr,:))==0) && (frameCurr>1))
                handles.normRolls(frameCurr,:) = handles.normRolls(frameCurr-1,:) ;
            
                %disp('Using roll vector from previous frame') ;
            end
            
            % if still zero, start with zero roll
            if (norm(handles.normRolls(frameCurr,:))==0)
                handles.normRolls(frameCurr,:) = - cross( rollHat, [0 0 1]) ;
                %disp('starting from rho=0') ;
            end
            normRoll = handles.normRolls(frameCurr,:) ;
            % make sure the roll vector is perpendicular to AHat  and renormalize
            normRoll = normRoll - rollHat * dot(normRoll, rollHat) ;
            normRoll = normRoll / norm(normRoll) ;
            handles.normRolls(frameCurr,:) = normRoll ;
            
            stubvec  = handles.Lstub * RotatePoint(normRoll,[0 0 0],rollHat,90) ;
            psiHat = handles.psiHats(frameCurr,:) ;
            
            %============================
            % get wing vector info
            %============================
            cr = handles.rightWingCM(frameCurr,:) ;
            cl = handles.leftWingCM(frameCurr,:) ;
            
            span1Hat = handles.span1Hats(frameCurr,:) ;
            span2Hat = handles.span2Hats(frameCurr,:) ;
            chord1Hat = handles.chord1Hats(frameCurr,:) ;
            chord2Hat = handles.chord2Hats(frameCurr,:) ;
            %--------------------------------------------------------------------------
            
            %============================
            % voxel plotting
            %============================
            if handles.bodyChangeFlag
                set(handles.hBody, 'XData', coords(bodyRows,1),...
                    'YData', coords(bodyRows,2),...
                    'ZData', coords(bodyRows,3));
            end
            if handles.rightWingChangeFlag
                set(handles.hRightWing, 'XData', coords(rightWingRows,1),...
                    'YData', coords(rightWingRows,2),...
                    'ZData', coords(rightWingRows,3));
            end
            if handles.leftWingChangeFlag
                set(handles.hLeftWing, 'XData', coords(leftWingRows,1),...
                    'YData', coords(leftWingRows,2),...
                    'ZData', coords(leftWingRows,3));
            end
            
            %============================
            % body vector plotting
            %============================
            if handles.AhatChangeFlag
                x_data = [cb(1), cb(1)+ scale*15*rollHat(1)] ;
                y_data = [cb(2), cb(2)+ scale*15*rollHat(2)] ;
                z_data = [cb(3), cb(3)+ scale*15*rollHat(3)] ;
            
                set(handles.hAhat, 'XData', x_data,...
                    'YData', y_data,...
                    'ZData', z_data);
            end
            if handles.stublineChangeFlag
            
                x_data = [cb(1), cb(1)+ stubvec(1)] ;
                y_data = [cb(2), cb(2)+ stubvec(2)] ;
                z_data = [cb(3), cb(3)+ stubvec(3)] ;
                set(handles.hStubline, 'XData', x_data,...
                    'YData', y_data,...
                    'ZData', z_data);
            end
            if handles.nrChangeFlag
                x_data = [cb(1)+stubvec(1)-Lbar*normRoll(1), cb(1)+stubvec(1)+Lbar*normRoll(1)] ;
                y_data = [cb(2)+stubvec(2)-Lbar*normRoll(2), cb(2)+stubvec(2)+Lbar*normRoll(2)] ;
                z_data = [cb(3)+stubvec(3)-Lbar*normRoll(3), cb(3)+stubvec(3)+Lbar*normRoll(3)] ;
                set(handles.hNr, 'XData', x_data,...
                    'YData', y_data,...
                    'ZData', z_data);
            end
            if handles.phChangeFlag
                x_data = [cb(1),cb(1)+scale*7*psiHat(1)] ;
                y_data = [cb(2),cb(2)+scale*7*psiHat(2)] ;
                z_data = [cb(3),cb(3)+scale*7*psiHat(3)] ;
                set(handles.hPh, 'XData', x_data,...
                    'YData', y_data,...
                    'ZData', z_data);
            end
            
            %============================
            % wing vector plotting
            %============================
            if handles.spanRightChangeFlag
                x_data = [cr(1), cr(1)+ scale*10*span1Hat(1)] ;
                y_data = [cr(2), cr(2)+ scale*10*span1Hat(2)] ;
                z_data = [cr(3), cr(3)+ scale*10*span1Hat(3)] ;
            
                set(handles.hSpanRight, 'XData', x_data,...
                    'YData', y_data,...
                    'ZData', z_data);
            end
            if handles.chordRightChangeFlag
                x_data = [cr(1), cr(1)+ scale*8*chord1Hat(1)] ;
                y_data = [cr(2), cr(2)+ scale*8*chord1Hat(2)] ;
                z_data = [cr(3), cr(3)+ scale*8*chord1Hat(3)] ;
                set(handles.hChordRight, 'XData', x_data,...
                    'YData', y_data,...
                    'ZData', z_data);
            end
            if handles.spanLeftChangeFlag
                x_data = [cl(1), cl(1)+ scale*10*span2Hat(1)] ;
                y_data = [cl(2), cl(2)+ scale*10*span2Hat(2)] ;
                z_data = [cl(3), cl(3)+ scale*10*span2Hat(3)] ;
                set(handles.hSpanLeft, 'XData', x_data,...
                    'YData', y_data,...
                    'ZData', z_data);
            end
            if handles.chordLeftChangeFlag
                x_data = [cl(1), cl(1)+ scale*8*chord2Hat(1)] ;
                y_data = [cl(2), cl(2)+ scale*8*chord2Hat(2)] ;
                z_data = [cl(3), cl(3)+ scale*8*chord2Hat(3)] ;
                set(handles.hChordLeft, 'XData', x_data,...
                    'YData', y_data,...
                    'ZData', z_data);
            end
            
            set(handles.frame_info,'String',...
                ['t(ms)=' num2str(app.tmsCurr) '  Frame # ' ...
                num2str(app.frameCurr) ' of ' num2str(handles.Nimages)])
            
            guidata(hObject, handles);
        end
        
        function resetInteractions(app, event)
            % This function resets the states of the toggle tools that
            % impact user interactions.  It also resets the figure interactions.
             
            % Find all tools to reset.  Exclude the tool associated
            % with the event.
            interactiveTools = [app.ui_zoom_in, app.ui_zoom_out, app.ui_rotate, app.ui_pan];
            interactiveTools(event.Source == interactiveTools) = [];
             
            % Set the state of the tools to 'off'.
            [interactiveTools.State] = deal('off');
             
            % Set figure interactions to 'off'.
            datacursormode(app.figure1, 'off')
            rotate3d(app.figure1, 'off');
            pan(app.figure1, 'off');
            zoom(app.figure1,'off');
        end

        %{
        function updateSliderWithButtons(app,handles,hObject,slider,adjustVal)
            frameCurr = handles.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;
            
            if newVal < slider.Limits(2) && newVal > slider.Limits(1)
                slider.Value = newVal;
                handles.bodyCM(frameCurr,1) = handles.bodyCM(frameCurr,1) + adjustVal;
                handles = updateChangeFlags(app,handles,'bodyvec');
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end
        %}
    end
    
    methods (Access = public)
        
        function EZViewReturn(app, tmsCurr, frameCurr, originalHandles)
            [hObject,~,~] = convertToGUIDECallbackArguments(app);
            %set(handles,originalHandles);
            
            app.tmsCurr = tmsCurr;
            app.frameCurr = frameCurr;

            handles = updateChangeFlags(app,originalHandles,'all');

            guidata(hObject,handles);

            resetToggleButtons(app, hObject, app.frameCurr) % reset buttons
            resetSliders(app, hObject, app.frameCurr)  % reset sliders

            guidata(hObject,handles);

            updateDisplay(app,hObject);

            guidata(hObject,handles);
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function correctionGUI_sam_OpeningFcn(app, varargin)
            % --- Executes just before correctionGUI_sam is made visible.
            
            % Ensure that the app appears on screen when run
            movegui(app.figure1, 'onscreen');
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app); %#ok<ASGLU>
            
            % This function has no output args, see OutputFcn.
            % hObject    handle to figure
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            % varargin   unrecognized PropertyName/PropertyValue pairs from the
            %            command line (see VARARGIN)
            tmsCurr_default = -10 ; %-10 ;
            saveEvery_default = 200000 ; %200
            %startdir_default = 'D:\Fly Data\VNC MN Chrimson\28_05042019\Analysis\Unsorted\Expr_28_mov_006\' ;%'D:\Fly Data\VNC Motor Lines\' ;
            startdir_default = 'Y:\Abby\Haltere_Experiment_Data\82_01042024\Analysis\Unsorted\Expr_82_mov_023' ;
            bodyFrameFlag_default = false ;
            largePertFlag_default = true ;
            
            % marker size value (different screens can make points harder to see)
            handles.marker_size = 4 ;  %3
            % assign default values if no inputs provided
            switch nargin
                case 1
                    app.tmsCurr = tmsCurr_default ;
                    handles.saveEvery = saveEvery_default ;
                    handles.startdir = startdir_default ;
                    handles.bodyFrameFlag = bodyFrameFlag_default ;
                    handles.largePertFlag = largePertFlag_default ;
            
                    handles.tmsCurr_default = tmsCurr_default ;
                case 2
                    app.tmsCurr = varargin{1} ;
                    handles.saveEvery = saveEvery_default ;
                    handles.startdir = startdir_default ;
                    handles.bodyFrameFlag = bodyFrameFlag_default ;
                    handles.largePertFlag = largePertFlag_default ;
            
                    handles.tmsCurr_default = varargin{1} ;
            
                case 3
                    app.tmsCurr = varargin{1} ;
                    handles.saveEvery = varargin{2} ;
                    handles.startdir = startdir_default ;
                    handles.bodyFrameFlag = bodyFrameFlag_default ;
                    handles.largePertFlag = largePertFlag_default ;
            
                    handles.tmsCurr_default = varargin{1} ;
                case 4
                    app.tmsCurr = varargin{1} ;
                    handles.saveEvery = varargin{2} ;
                    handles.startdir = varargin{3} ;
                    handles.bodyFrameFlag = bodyFrameFlag_default ;
                    handles.largePertFlag = largePertFlag_default ;
            
                    handles.tmsCurr_default = varargin{1} ;
                case 5
                    app.tmsCurr = varargin{1} ;
                    handles.saveEvery = varargin{2} ;
                    handles.startdir = varargin{3} ;
                    handles.bodyFrameFlag = varargin{4} ;
                    handles.largePertFlag = largePertFlag_default ;
            
                    handles.tmsCurr_default = varargin{1} ;
                case 6
                    app.tmsCurr = varargin{1} ;
                    handles.saveEvery = varargin{2} ;
                    handles.startdir = varargin{3} ;
                    handles.bodyFrameFlag = varargin{4} ;
                    handles.largePertFlag = varargin{6} ;
            
                    handles.tmsCurr_default = varargin{1} ;
                otherwise
                    disp('Wrong number of input arguments')
                    keyboard
            end
            % Choose default command line output for correctionGUI_sam
            handles.output = hObject;
            
            % initialize plot handles
            handles = initPlots(app, handles) ;
            
            % try to deal with zoom/keypress conflict
            
            handles.hManager = uigetmodemanager(handles.figure1);
            try
                set(handles.hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
            catch
                [handles.hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
            end
            set(handles.figure1, 'WindowKeyPressFcn',@figure1_WindowKeyPressFcn);
            set(handles.figure1, 'KeyPressFcn', []);
            
            
            % flag to let us know whether or not data is loaded
            handles.dataFlag = false ;
            handles.datapath_curr = [] ;
            handles.dirlist = [] ;
            
            % setting to determine whether or not we should recalculate angles every
            % time we save manually corrected data
            handles.calcNewAnglesFlag = true ;
            
            % reference structures to make callbacks easier
            adjustTagStruct = struct() ;
            %slider_handles = findall(0,'Style','Slider') ;
            slider_handles = [handles.body_roll, handles.body_pitch, handles.body_phi, handles.l_eta, handles.l_theta, handles.l_phi, handles.l_z,...
                handles.l_y, handles.l_x, handles.r_eta, handles.r_theta, handles.r_phi, handles.r_z, handles.rx, handles.body_z,...
                handles.body_y, handles.body_x, handles.r_y];
            for i = 1:length(slider_handles)
                adjustTagStruct.(slider_handles(i).Tag) = i + 1 ;
            end
            handles.adjustTagStruct = adjustTagStruct ;
            
            % Update handles structure
            guidata(hObject, handles);
        end

        % Button pushed function: back
        function back_Callback(app, event)
            % --- Executes on button press in back.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to back (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % shift current timing/frame
            app.frameCurr = max([app.frameCurr - 1, 1]) ;
            app.tmsCurr = handles.tvec(app.frameCurr) ;
            handles = updateChangeFlags(app, handles,'all') ;
            resetToggleButtons(app, hObject, app.frameCurr) % reset buttons
            resetSliders(app, hObject, app.frameCurr)  % reset sliders
            if mod(app.frameCurr,handles.saveEvery) == 0
                saveData(app, hObject)
            end
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Button pushed function: bback
        function bback_Callback(app, event)
            % --- Executes on button press in bback.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to bback (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % shift current timing/frame
            app.frameCurr = app.frameCurr - 10 ;
            app.tmsCurr = handles.tvec(app.frameCurr) ;
            handles = updateChangeFlags(app, handles,'all') ;
            resetToggleButtons(app, hObject, app.frameCurr)
            resetSliders(app, hObject, app.frameCurr)
            if mod(app.frameCurr,handles.saveEvery) == 0
                saveData(app, hObject)
            end
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: body_phi
        function body_phi_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to body_phi (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.rollHats(frameCurr,:) = ...
                RotatePoint(handles.rollHats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
            handles.normRolls(frameCurr,:) = ...
                RotatePoint(handles.normRolls(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
            handles.psiHats(frameCurr,:) = ...
                RotatePoint(handles.psiHats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
            handles = updateChangeFlags(app, handles,'bodyvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: body_pitch
        function body_pitch_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to body_pitch (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.rollHats(frameCurr,:) = ...
                RotatePoint2(handles.rollHats(frameCurr,:),[0 0 0],...
                handles.psiHats(frameCurr,:),adjust_val);
            handles.normRolls(frameCurr,:) = ...
                RotatePoint2(handles.normRolls(frameCurr,:),[0 0 0],...
                handles.psiHats(frameCurr,:),adjust_val);
            handles = updateChangeFlags(app, handles,'bodyvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: body_roll
        function body_roll_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to body_roll (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.normRolls(frameCurr,:) = ...
                RotatePoint(handles.normRolls(frameCurr,:),[0 0 0],...
                handles.rollHats(frameCurr,:),-1*adjust_val);
            handles = updateChangeFlags(app, handles,'bodyvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: body_x
        function body_x_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to body_x (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.bodyCM(frameCurr,1) = handles.bodyCM(frameCurr,1) + adjust_val ;
            handles = updateChangeFlags(app, handles,'bodyvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: body_y
        function body_y_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to body_x (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.bodyCM(frameCurr,2) = handles.bodyCM(frameCurr,2) + adjust_val ;
            handles = updateChangeFlags(app, handles,'bodyvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: body_z
        function body_z_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to body_z (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.bodyCM(frameCurr,3) = handles.bodyCM(frameCurr,3) + adjust_val ;
            handles = updateChangeFlags(app, handles,'bodyvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Button pushed function: clear_data
        function clear_data_Callback(app, event)
            % --- Executes on button press in clear_data.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to clear_data (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            tempdirlist = handles.dirlist ;
            remove_idx = arrayfun(@(x) strcmp(x.name, handles.datapath_curr), ...
                tempdirlist) ;
            
            % take data name out of list box (if we have things to remove)
            if (sum(remove_idx) > 0)
                set(handles.data_dir,'String',{tempdirlist(~remove_idx).name})
                % take data file out of directory
                handles.dirlist = tempdirlist(~remove_idx) ;
            end
            
            clearDisplay(app, hObject) ;  % reset plot window
            handles.dataFlag = false ;      % signal that there's no more data
            handles.datapath_curr = [] ;    % empty current data path
            
            guidata(hObject, handles);
        end

        % Button pushed function: cluster_wings
        function cluster_wings_Callback(app, event)
            % --- Executes on button press in cluster_wings.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to cluster_wings (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            % -----------------------------------
            % read in initial data
            frameCurr = app.frameCurr ;
            N_vox_min = 200 ;
            
            row1 = handles.frameStartInd(frameCurr) ;
            row2 = handles.frameEndInd(frameCurr) ;
            IDX = handles.RESIDX(row1:row2,:) ;  % (isbody, isrightwing, isleftwing)
            rightWingRows = (IDX(:,handles.rightWingInd)==1) ;
            leftWingRows  = (IDX(:,handles.leftWingInd)==1) ;
            
            N_vox_R = sum(rightWingRows) ; % number of voxels in each wing
            N_vox_L = sum(leftWingRows) ;
            
            % --------------------------------------------------------
            % smooth/interpolate wing CMs and tips
            %
            [rightWingCM_interp, ~, rightWingTips_interp] = ...
                interpolateWingCM(handles,'right') ;
            [leftWingCM_interp, ~, leftWingTips_interp] = ...
                interpolateWingCM(handles,'left') ;
            % --------------------------------------------------------
            % cluster the voxels for whichever wing has more voxels
            if (N_vox_R > N_vox_L) && (N_vox_R > N_vox_min)
                wing_str = 'right' ;
            elseif (N_vox_R < N_vox_L) && (N_vox_L > N_vox_min)
                wing_str = 'left' ;
            else
                disp('Uncertain which wing to cluster')
                return
            end
            % -------------------------------------------------------------------
            % perform clustering
            [wingVox, wingRows, label_idx, centroids, badClusterFlag] = ...
                clusterWings(handles, frameCurr, wing_str) ;
            
            % ------------------------------------------------------------------
            % did clustering work?
            if ~badClusterFlag
                % if we did successfully cluster, we now need to figure out
                % which is left vs right. we'll do this by checking against the
                % interpolated center of mass
                right_centroid_dist = myNorm(centroids -  ...
                    repmat(rightWingCM_interp(frameCurr,:),2,1)) ;
                left_centroid_dist = myNorm(centroids -  ...
                    repmat(leftWingCM_interp(frameCurr,:),2,1)) ;
                [~, right_idx] = min(right_centroid_dist) ;
                [~, left_idx] = min(left_centroid_dist) ;
            
                % in case both blobs think they should be the same wing
                if (right_idx == left_idx)
                    disp('could not determine left vs right cluster')
                else
                    % otherwise, assign the new wing data to our output struct
                    wingVoxR = wingVox(label_idx == right_idx,:) ;
                    rightCM = centroids(right_idx,:) ;
                    wingRowsR = wingRows(:,right_idx) ;
                    wingVoxL = wingVox(label_idx == left_idx,:) ;
                    leftCM = centroids(left_idx,:) ;
                    wingRowsL = wingRows(:,left_idx) ;
            
                    % now calculate wing vectors
                    rightRefVecs = [handles.bodyCM(frameCurr,:); ...
                        handles.bodyCM(frameCurr-1,:); ...
                        rightWingTips_interp(frameCurr-1,:)] ;
                    [spanHatR, chordHatR, chordAltHatR, ~, wingTipR] = ...
                        estimate_wing_vecs(wingVoxR, rightRefVecs, handles.wingLength, ...
                        [], [], rightCM) ;
            
                    leftRefVecs = [handles.bodyCM(frameCurr,:) ; ...
                        handles.bodyCM(frameCurr-1,:); ...
                        leftWingTips_interp(frameCurr-1,:)] ;
                    [spanHatL, chordHatL, chordAltHatL, ~, wingTipL] = ...
                        estimate_wing_vecs(wingVoxL, leftRefVecs, handles.wingLength, ...
                        [], [], leftCM) ;
            
                    % add to storage arrays
                    handles.rightWingCM(frameCurr,:) = rightCM ;
                    handles.span1Hats(frameCurr,:) = spanHatR ;
                    handles.chord1Hats(frameCurr,:) = chordHatR ;
                    %handles.rightChordAltHats(frameCurr,:) = chordAltHatR ;
                    handles.rightWingTips(frameCurr,:) = wingTipR ;
            
                    handles.leftWingCM(frameCurr,:) = leftCM ;
                    handles.span2Hats(frameCurr,:) = spanHatL ;
                    handles.chord2Hats(frameCurr,:) = chordHatL ;
                    %handles.leftChordAltHats(frameCurr,:) = chordAltHatL ;
                    handles.leftWingTips(frameCurr,:) = wingTipL ;
            
                    %...including fucking voxels, blurg
                    handles.RESIDX(row1:row2,2:3) = [wingRowsR, wingRowsL] ;
            
                    % now update display and save changes to data
                    handles = updateChangeFlags(app, handles,'all') ;
                    guidata(hObject, handles);
                    updateDisplay(app, hObject);
                end
            else
                disp('clustering failed')
            end
        end

        % Value changed function: data_dir
        function data_dir_Callback(app, event)
            % --- Executes on selection change in data_dir.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to data_dir (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: contents = cellstr(get(hObject,'String')) returns data_dir contents as cell array
            %        contents{get(hObject,'Value')} returns selected item from data_dir
            index_selected = get(hObject,'Value');
            contents = get(hObject,'String');
            item_selected = contents{index_selected};
            
            handles.datapath_curr = item_selected ;
            
            guidata(hObject, handles);
        end

        % Button pushed function: ffwd
        function ffwd_Callback(app, event)
            % --- Executes on button press in ffwd.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to ffwd (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            if (app.frameCurr + 10) >= handles.Nimages
                saveData(app, hObject)
                disp('Not enough remaining frames to skip forward')
                return
            end
            
            % shift current timing/frame
            app.frameCurr = app.frameCurr + 10 ;
            app.tmsCurr = handles.tvec(app.frameCurr) ;
            handles = updateChangeFlags(app, handles,'all') ;
            resetToggleButtons(app, hObject, app.frameCurr) % reset buttons
            resetSliders(app, hObject, app.frameCurr)  % reset sliders
            if mod(app.frameCurr,handles.saveEvery) == 0
                saveData(app, hObject)
            end
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Window key press function: figure1
        function figure1_WindowKeyPressFcn(app, event)
            %beep
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to figure1 (see GCBO)
            % eventdata  structure with the following fields (see FIGURE)
            %	Key: name of the key that was pressed, in lower case
            %	Character: character interpretation of the key(s) that was pressed
            %	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
            % handles    structure with handles and user data (see GUIDATA)
            handles = guidata(hObject);
            disp(eventdata.Key)
            if handles.dataFlag
                switch eventdata.Key
                    case 'v'
                        uicontrol(handles.bback);
                        bback_Callback(app, event);
                    case 'b'
                        uicontrol(handles.back);
                        back_Callback(app, event);
                    case 'return'
                        uicontrol(handles.fwd);
                        fwd_Callback(app, event);
                    case 'f'
                        uicontrol(handles.ffwd);
                        ffwd_Callback(app, event);
                    case 'downarrow'
                        uicontrol(handles.bback);
                        bback_Callback(app,event);
                    case 'leftarrow'
                        uicontrol(handles.back);
                        back_Callback(app,event);
                    case 'rightarrow'
                        uicontrol(handles.fwd);
                        fwd_Callback(app,event);
                    case 'uparrow'
                        uicontrol(handles.ffwd);
                        ffwd_Callback(app,event);
                    case 's'
                        uicontrol(handles.bback);
                        bback_Callback(app,event);
                    case 'a'
                        uicontrol(handles.back);
                        back_Callback(app,event);
                    case 'd'
                        uicontrol(handles.fwd);
                        fwd_Callback(app,event);
                    case 'w'
                        uicontrol(handles.ffwd);
                        ffwd_Callback(app,event);
                    case '1'
                        %uicontrol(handles.phi_view)
                        phi_viewButtonPushed(app,event)
                    case '2'
                        %uicontrol(handles.r_theta_view)
                        r_theta_viewButtonPushed(app,event)
                    case '3'
                        %uicontrol(handles.l_theta_view)
                        l_theta_viewButtonPushed(app,event)
                    case '4'
                        %uicontrol(handles.r_eta_view)
                        r_eta_viewButtonPushed(app,event)
                    case '5'
                        %uicontrol(handles.l_eta_view)
                        l_eta_viewButtonPushed(app,event)
                    case '6'
                        %uicontrol(handles.roll_view)
                        roll_viewButtonPushed(app,event)
                end
            end
            
        end

        % Button pushed function: fwd
        function fwd_Callback(app, event)
            % --- Executes on button press in fwd.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            % hObject    handle to fwd (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            if app.frameCurr >= handles.Nimages
                saveData(app, hObject)
                disp('Movie completed!')
                return
            end
            
            % shift current timing/frame
            app.frameCurr = app.frameCurr + 1 ;
            app.tmsCurr = handles.tvec(app.frameCurr) ;
            handles = updateChangeFlags(app, handles,'all') ;
            resetToggleButtons(app, hObject, app.frameCurr) % reset buttons
            resetSliders(app, hObject, app.frameCurr)  % reset sliders
            if mod(app.frameCurr,handles.saveEvery) == 0
                saveData(app, hObject)
            end
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: ignoreFrame
        function ignoreFrame_Callback(app, event)
            % --- Executes on button press in ignoreFrame.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to ignoreFrame (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            button_state = get(hObject,'Value');
            frameCurr = app.frameCurr ;
            
            if button_state == get(hObject,'Max')
                handles.ignoreFlag(frameCurr) = true ;
                handles.ignoreFrames = sort(unique([handles.ignoreFrames, frameCurr])) ;
                set(hObject,'BackgroundColor','red')
            elseif button_state == get(hObject,'Min')
                handles.ignoreFlag(frameCurr) = false ;
                handles.ignoreFrames = ...
                    sort(unique(handles.ignoreFrames(handles.ignoreFrames ~= frameCurr))) ;
                defaultColor = get(0,'defaultUicontrolBackgroundColor') ;
                set(hObject,'BackgroundColor',defaultColor)
            end
            guidata(hObject, handles);
        end

        % Value changed function: l_eta
        function l_eta_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to l_eta (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.chord2Hats(frameCurr,:) = ...
                RotatePoint(handles.chord2Hats(frameCurr,:),[0 0 0],...
                handles.span2Hats(frameCurr,:),-1*adjust_val);
            handles = updateChangeFlags(app, handles,'lchord') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: l_phi
        function l_phi_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to l_phi (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            
            handles.span2Hats(frameCurr,:) = ...
                RotatePoint(handles.span2Hats(frameCurr,:),[0 0 0],[0 0 1],-1*adjust_val);
            handles.chord2Hats(frameCurr,:) = ...
                RotatePoint(handles.chord2Hats(frameCurr,:),[0 0 0],[0 0 1],-1*adjust_val);
            %chord2Hat2s(i,:) = RotatePoint(chord2Hat2s(i,:),[0 0 0],[0 0 1],-10);
            handles.phi2Hats(frameCurr,:) = ...
                RotatePoint(handles.phi2Hats(frameCurr,:),[0 0 0],[0 0 1],-1*adjust_val);
            handles = updateChangeFlags(app, handles,'lvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: l_theta
        function l_theta_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to l_theta (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.span2Hats(frameCurr,:) = ...
                RotatePoint2(handles.span2Hats(frameCurr,:),[0 0 0],...
                handles.phi2Hats(frameCurr,:),-1*adjust_val);
            handles.chord2Hats(frameCurr,:) = ...
                RotatePoint2(handles.chord2Hats(frameCurr,:),[0 0 0],...
                handles.phi2Hats(frameCurr,:),-1*adjust_val);
            %chord2Hat2s(i,:) = RotatePoint2(chord2Hat2s(i,:),[0 0 0],phi2Hats(i,:),10);
            handles = updateChangeFlags(app, handles,'lvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: l_x
        function l_x_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to l_x (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.leftWingCM(frameCurr,1) = handles.leftWingCM(frameCurr,1) + adjust_val ;
            handles = updateChangeFlags(app, handles,'lvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: l_y
        function l_y_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to l_y (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.leftWingCM(frameCurr,2) = handles.leftWingCM(frameCurr,2) + adjust_val ;
            handles = updateChangeFlags(app, handles,'lvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: l_z
        function l_z_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to l_z (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.leftWingCM(frameCurr,3) = handles.leftWingCM(frameCurr,3) + adjust_val ;
            handles = updateChangeFlags(app, handles,'lvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Button pushed function: load_data
        function load_data_Callback(app, event)
            % --- Executes on button press in load_data.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to load_data (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            if handles.dataFlag
                clearDisplay(app, hObject) ;
            end
            handles = loadData(app, hObject, handles) ;
        end

        % Button pushed function: open_data_dir
        function open_data_dir_Callback(app, event)
            % --- Executes on button press in open_data_dir.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to open_data_dir (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % -----------------------------------------------
            % see if there's a set directory to start search
            try
                startdir = handles.startdir;
            catch
                startdir = pwd;
            end
            % ----------------------------------------------------------------------
            % if there's not already data in the directory, just automatically load
            if isfield(handles, 'dirlist')
                if ~isempty(handles.dirlist)
                    loadFlag = false ;
                else
                    loadFlag = true ;
                end
            else
                loadFlag = true ;
            end
            % ---------------------------------------------
            % use file picker to get files
            tempdirlist = uipickfiles('filter',startdir, 'output', 'struct');
            if ~isempty(tempdirlist)
                mat_file_ind = arrayfun(@(x) contains(x.name,'.mat') & ...
                    contains(x.name,'Expr'),tempdirlist) ;
                tempdirlist = tempdirlist(mat_file_ind) ;
                tempdirlist = vertcat(handles.dirlist, tempdirlist) ;
                handles.dirlist = tempdirlist;
                set(handles.data_dir,'String',{tempdirlist(~[tempdirlist(:).isdir]).name})
                handles.datapath_curr = tempdirlist(1).name ;
            end
            if loadFlag && ~isempty(handles.datapath_curr)
                handles = loadData(app, hObject, handles) ;
            else
                guidata(hObject, handles);
            end
        end

        % Value changed function: r_eta
        function r_eta_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to r_eta (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.chord1Hats(frameCurr,:) = ...
                RotatePoint(handles.chord1Hats(frameCurr,:),[0 0 0],...
                handles.span1Hats(frameCurr,:),-1*adjust_val);
            handles = updateChangeFlags(app, handles,'rchord') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: r_phi
        function r_phi_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to r_phi (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            
            handles.span1Hats(frameCurr,:) = ...
                RotatePoint(handles.span1Hats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
            handles.chord1Hats(frameCurr,:) = ...
                RotatePoint(handles.chord1Hats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
            %chord1Hat2s(i,:) = RotatePoint(chord1Hat2s(i,:),[0 0 0],[0 0 1],-10);
            handles.phi1Hats(frameCurr,:) = ...
                RotatePoint(handles.phi1Hats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
            handles = updateChangeFlags(app, handles,'rvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: r_theta
        function r_theta_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to r_theta (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.span1Hats(frameCurr,:) = ...
                RotatePoint2(handles.span1Hats(frameCurr,:),[0 0 0],...
                handles.phi1Hats(frameCurr,:),-1*adjust_val);
            handles.chord1Hats(frameCurr,:) = ...
                RotatePoint2(handles.chord1Hats(frameCurr,:),[0 0 0],...
                handles.phi1Hats(frameCurr,:),-1*adjust_val);
            %chord1Hat2s(i,:) = RotatePoint2(chord1Hat2s(i,:),[0 0 0],phi1Hats(i,:),10);
            handles = updateChangeFlags(app, handles,'rvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: r_y
        function r_y_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to r_y (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.rightWingCM(frameCurr,2) = handles.rightWingCM(frameCurr,2) + adjust_val ;
            handles = updateChangeFlags(app, handles,'rvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: r_z
        function r_z_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to r_z (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.rightWingCM(frameCurr,3) = handles.rightWingCM(frameCurr,3) + adjust_val ;
            handles = updateChangeFlags(app, handles,'rvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Value changed function: rx
        function rx_Callback(app, event)
            % --- Executes on slider movement.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to rx (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            frameCurr = app.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.rightWingCM(frameCurr,1) = handles.rightWingCM(frameCurr,1) + adjust_val ;
            handles = updateChangeFlags(app, handles,'rvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Button pushed function: save_manualCorrRangeMS
        function save_manualCorrRangeMS_Callback(app, event)
            % --- Executes on button press in save_manualCorrRangeMS.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to save_manualCorrRangeMS (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            prompt = {'Enter manual correction start time (ms):',...
                'Enter manual correction end time (ms):'};
            dlgtitle = 'Save Manual correction range';
            dims = [1 40];
            definput = {'-10','30'};
            userInput = inputdlg(prompt,dlgtitle,dims,definput) ;
            
            if ~isempty(userInput)
                manualCorrRangeMS = cellfun(@(y) str2double(y), userInput) ;
                if size(manualCorrRangeMS,1) > size(manualCorrRangeMS,2)
                    manualCorrRangeMS = manualCorrRangeMS' ;
                end
                handles.manualCorrRangeMS = manualCorrRangeMS ;
            else
                handles.manualCorrRangeMS = [] ;
            end
            guidata(hObject, handles);
        end

        % Value changed function: save_roll
        function save_roll_Callback(app, event)
            % --- Executes on button press in save_roll.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to save_roll (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            button_state = get(hObject,'Value');
            frameCurr = app.frameCurr ;
            
            if button_state == get(hObject,'Max')
                handles.rhoFlag(frameCurr) = true ;
                handles.rhoTimes = sort(unique([handles.rhoTimes, frameCurr])) ;
                set(hObject,'BackgroundColor','green')
            elseif button_state == get(hObject,'Min')
                handles.rhoFlag(frameCurr) = false ;
                handles.rhoTimes = ...
                    sort(unique(handles.rhoTimes(handles.rhoTimes ~= frameCurr))) ;
                defaultColor = get(0,'defaultUicontrolBackgroundColor') ;
                set(hObject,'BackgroundColor',defaultColor)
            end
            guidata(hObject, handles);
        end

        % Button pushed function: swap_wings
        function swap_wings_Callback(app, event)
            % --- Executes on button press in swap_wings.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to swap_wings (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            frameCurr = app.frameCurr ;
            
            % swap voxels
            row1 = handles.frameStartInd(frameCurr) ;
            row2 = handles.frameEndInd(frameCurr) ;
            IDX = handles.RESIDX(row1:row2,:) ;  % (isbody, isrightwing, isleftwing)
            rightWingRows = (IDX(:,handles.rightWingInd)==1) ;
            leftWingRows  = (IDX(:,handles.leftWingInd)==1) ;
            allRows = row1:row2 ;
            
            handles.RESIDX(allRows(rightWingRows),handles.rightWingInd) = 0 ;
            handles.RESIDX(allRows(rightWingRows),handles.leftWingInd) = 1 ;
            handles.RESIDX(allRows(leftWingRows),handles.rightWingInd) = 1 ;
            handles.RESIDX(allRows(leftWingRows),handles.leftWingInd) = 0 ;
            
            % swap vectors
            span1Hats = handles.span1Hats ;
            span2Hats = handles.span2Hats ;
            chord1Hats = handles.chord1Hats ;
            chord2Hats = handles.chord2Hats ;
            phi1Hats = handles.phi1Hats ;
            phi2Hats = handles.phi2Hats ;
            
            handles.span1Hats(frameCurr,:) = span2Hats(frameCurr,:) ;
            handles.chord1Hats(frameCurr,:) = chord2Hats(frameCurr,:) ;
            handles.phi1Hats(frameCurr,:) = phi2Hats(frameCurr,:) ;
            
            handles.span2Hats(frameCurr,:) = span1Hats(frameCurr,:) ;
            handles.chord2Hats(frameCurr,:) = chord1Hats(frameCurr,:) ;
            handles.phi2Hats(frameCurr,:) = phi1Hats(frameCurr,:) ;
            
            % swap cm
            rightWingCM = handles.rightWingCM ;
            leftWingCM = handles.leftWingCM ;
            
            handles.rightWingCM(frameCurr,:) = leftWingCM(frameCurr,:) ;
            handles.leftWingCM(frameCurr,:) = rightWingCM(frameCurr,:) ;
            
            % update structure and plot
            handles.wingSwapFlag = ~handles.wingSwapFlag ;
            handles = updateChangeFlags(app, handles,'all') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);
        end

        % Callback function: ui_save
        function ui_save_ClickedCallback(app, event)
            % --------------------------------------------------------------------
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to ui_save (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % --- Executes during object creation, after setting all properties.
            saveData(app, hObject)
        end

        % Callback function: ui_pan
        function ui_pan_ClickedCallback(app, event)
            % Reset the states of interactive tools and reset all figure
            % interactions.
            app.resetInteractions(event);
             
            % Enable or disable pan based on the
            % tool's current state.
            state = app.ui_pan.State;
            pan(app.figure1, char(state));
        end

        % Callback function: ui_rotate
        function ui_rotate_ClickedCallback(app, event)
            % Reset the states of interactive tools and reset all figure
            % interactions.
            app.resetInteractions(event);
             
            % Enable or disable rotation based on the
            % tool's current state.
            state = app.ui_rotate.State;
            rotate3d(app.figure1, char(state));
        end

        % Callback function: ui_zoom_out
        function ui_zoom_out_ClickedCallback(app, event)
            % Reset the states of interactive tools and reset all figure
            % interactions.
            app.resetInteractions(event);
             
            % Enable or disable zoom-out based on the
            % tool's current state.
            state = app.ui_zoom_out.State;
            zoomModeObject = zoom(app.figure1);
            if state
                zoomModeObject.Direction = 'out';
                zoomModeObject.Enable = 'on';
            else
                zoomModeObject.Enable = 'off';
            end
        end

        % Callback function: ui_zoom_in
        function ui_zoom_in_ClickedCallback(app, event)
            % Reset the states of interactive tools and reset all figure
            % interactions.
            app.resetInteractions(event);
             
            % Enable or disable zoom-in based on the
            % tool's current state.
            state = app.ui_zoom_in.State;
            zoomModeObject = zoom(app.figure1);
            if state
                zoomModeObject.Direction = 'in';
                zoomModeObject.Enable = 'on';
            else
                zoomModeObject.Enable = 'off';
            end
        end

        % Button pushed function: body_x_inc
        function body_x_incButtonPushed(app, event)
            %{
            frameCurr = handles.frameCurr ;
            [handles, adjust_val] = sliderAdjust(app, hObject, handles) ;
            handles.bodyCM(frameCurr,1) = handles.bodyCM(frameCurr,1) + adjust_val ;
            handles = updateChangeFlags(app, handles,'bodyvec') ;
            guidata(hObject, handles);
            updateDisplay(app, hObject);

            %}
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_x;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;
            
            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.bodyCM(frameCurr,1) = handles.bodyCM(frameCurr,1) + adjustVal;
                handles = updateChangeFlags(app,handles,'bodyvec');
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_x_dec
        function body_x_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_x;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;
            
            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.bodyCM(frameCurr,1) = handles.bodyCM(frameCurr,1) + adjustVal;
                handles = updateChangeFlags(app,handles,'bodyvec');
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_y_dec
        function body_y_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_y;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;
            
            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.bodyCM(frameCurr,2) = handles.bodyCM(frameCurr,2) + adjustVal;
                handles = updateChangeFlags(app,handles,'bodyvec');
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_y_inc
        function body_y_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_y;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;
            
            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.bodyCM(frameCurr,2) = handles.bodyCM(frameCurr,2) + adjustVal;
                handles = updateChangeFlags(app,handles,'bodyvec');
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_z_dec
        function body_z_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_z;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;
            
            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.bodyCM(frameCurr,3) = handles.bodyCM(frameCurr,3) + adjustVal;
                handles = updateChangeFlags(app,handles,'bodyvec');
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_phi_dec
        function body_phi_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_phi;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rollHats(frameCurr,:) = ...
                    RotatePoint(handles.rollHats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles.normRolls(frameCurr,:) = ...
                    RotatePoint(handles.normRolls(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles.psiHats(frameCurr,:) = ...
                    RotatePoint(handles.psiHats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles = updateChangeFlags(app,handles,'bodyvec');
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end

        end

        % Button pushed function: body_pitch_dec
        function body_pitch_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_pitch;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rollHats(frameCurr,:) = ...
                    RotatePoint2(handles.rollHats(frameCurr,:),[0 0 0],...
                    handles.psiHats(frameCurr,:),-adjustVal);
                handles.normRolls(frameCurr,:) = ...
                    RotatePoint2(handles.normRolls(frameCurr,:),[0 0 0],...
                    handles.psiHats(frameCurr,:),-adjustVal);
                handles = updateChangeFlags(app, handles,'bodyvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_roll_dec
        function body_roll_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_roll;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.normRolls(frameCurr,:) = ...
                    RotatePoint(handles.normRolls(frameCurr,:),[0 0 0],...
                    handles.rollHats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'bodyvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: rx_dec
        function rx_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.rx;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rightWingCM(frameCurr,1) = handles.rightWingCM(frameCurr,1) + adjustVal ;
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_y_dec
        function r_y_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_y;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rightWingCM(frameCurr,2) = handles.rightWingCM(frameCurr,2) + adjustVal ;
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_z_dec
        function r_z_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_z;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rightWingCM(frameCurr,3) = handles.rightWingCM(frameCurr,3) + adjustVal ;
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_phi_dec
        function r_phi_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_phi;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.span1Hats(frameCurr,:) = ...
                    RotatePoint(handles.span1Hats(frameCurr,:),[0 0 0],[0 0 1],adjustVal);
                handles.chord1Hats(frameCurr,:) = ...
                    RotatePoint(handles.chord1Hats(frameCurr,:),[0 0 0],[0 0 1],adjustVal);
                handles.phi1Hats(frameCurr,:) = ...
                    RotatePoint(handles.phi1Hats(frameCurr,:),[0 0 0],[0 0 1],adjustVal);
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_theta_dec
        function r_theta_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_theta;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.span1Hats(frameCurr,:) = ...
                    RotatePoint2(handles.span1Hats(frameCurr,:),[0 0 0],...
                    handles.phi1Hats(frameCurr,:),-1*adjustVal);
                handles.chord1Hats(frameCurr,:) = ...
                    RotatePoint2(handles.chord1Hats(frameCurr,:),[0 0 0],...
                    handles.phi1Hats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_eta_dec
        function r_eta_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_eta;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.chord1Hats(frameCurr,:) = ...
                    RotatePoint(handles.chord1Hats(frameCurr,:),[0 0 0],...
                    handles.span1Hats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'rchord') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_x_dec
        function l_x_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_x;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.leftWingCM(frameCurr,1) = handles.leftWingCM(frameCurr,1) + adjustVal ;
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_y_dec
        function l_y_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_y;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.leftWingCM(frameCurr,2) = handles.leftWingCM(frameCurr,2) + adjustVal ;
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_z_dec
        function l_z_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_z;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.leftWingCM(frameCurr,3) = handles.leftWingCM(frameCurr,3) + adjustVal ;
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_phi_dec
        function l_phi_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_phi;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.span2Hats(frameCurr,:) = ...
                    RotatePoint(handles.span2Hats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles.chord2Hats(frameCurr,:) = ...
                    RotatePoint(handles.chord2Hats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles.phi2Hats(frameCurr,:) = ...
                    RotatePoint(handles.phi2Hats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_theta_dec
        function l_theta_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_theta;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.span2Hats(frameCurr,:) = ...
                    RotatePoint2(handles.span2Hats(frameCurr,:),[0 0 0],...
                    handles.phi2Hats(frameCurr,:),-1*adjustVal);
                handles.chord2Hats(frameCurr,:) = ...
                    RotatePoint2(handles.chord2Hats(frameCurr,:),[0 0 0],...
                    handles.phi2Hats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_eta_dec
        function l_eta_decButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_eta;
            adjustVal = -1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.chord2Hats(frameCurr,:) = ...
                    RotatePoint(handles.chord2Hats(frameCurr,:),[0 0 0],...
                    handles.span2Hats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'lchord') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_z_inc
        function body_z_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_z;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;
            
            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.bodyCM(frameCurr,3) = handles.bodyCM(frameCurr,3) + adjustVal;
                handles = updateChangeFlags(app,handles,'bodyvec');
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_phi_inc
        function body_phi_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_phi;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rollHats(frameCurr,:) = ...
                    RotatePoint(handles.rollHats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles.normRolls(frameCurr,:) = ...
                    RotatePoint(handles.normRolls(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles.psiHats(frameCurr,:) = ...
                    RotatePoint(handles.psiHats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles = updateChangeFlags(app,handles,'bodyvec');
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_pitch_inc
        function body_pitch_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_pitch;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rollHats(frameCurr,:) = ...
                    RotatePoint2(handles.rollHats(frameCurr,:),[0 0 0],...
                    handles.psiHats(frameCurr,:),-adjustVal);
                handles.normRolls(frameCurr,:) = ...
                    RotatePoint2(handles.normRolls(frameCurr,:),[0 0 0],...
                    handles.psiHats(frameCurr,:),-adjustVal);
                handles = updateChangeFlags(app, handles,'bodyvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: body_roll_inc
        function body_roll_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.body_roll;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.normRolls(frameCurr,:) = ...
                    RotatePoint(handles.normRolls(frameCurr,:),[0 0 0],...
                    handles.rollHats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'bodyvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: rx_inc
        function rx_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.rx;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rightWingCM(frameCurr,1) = handles.rightWingCM(frameCurr,1) + adjustVal ;
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_y_inc
        function r_y_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_y;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rightWingCM(frameCurr,2) = handles.rightWingCM(frameCurr,2) + adjustVal ;
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_z_inc
        function r_z_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_z;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.rightWingCM(frameCurr,3) = handles.rightWingCM(frameCurr,3) + adjustVal ;
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_phi_inc
        function r_phi_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_phi;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.span1Hats(frameCurr,:) = ...
                    RotatePoint(handles.span1Hats(frameCurr,:),[0 0 0],[0 0 1],adjustVal);
                handles.chord1Hats(frameCurr,:) = ...
                    RotatePoint(handles.chord1Hats(frameCurr,:),[0 0 0],[0 0 1],adjustVal);
                handles.phi1Hats(frameCurr,:) = ...
                    RotatePoint(handles.phi1Hats(frameCurr,:),[0 0 0],[0 0 1],adjustVal);
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_theta_inc
        function r_theta_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_theta;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.span1Hats(frameCurr,:) = ...
                    RotatePoint2(handles.span1Hats(frameCurr,:),[0 0 0],...
                    handles.phi1Hats(frameCurr,:),-1*adjustVal);
                handles.chord1Hats(frameCurr,:) = ...
                    RotatePoint2(handles.chord1Hats(frameCurr,:),[0 0 0],...
                    handles.phi1Hats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'rvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: r_eta_inc
        function r_eta_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.r_eta;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.chord1Hats(frameCurr,:) = ...
                    RotatePoint(handles.chord1Hats(frameCurr,:),[0 0 0],...
                    handles.span1Hats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'rchord') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_x_inc
        function l_x_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_x;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.leftWingCM(frameCurr,1) = handles.leftWingCM(frameCurr,1) + adjustVal ;
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_y_inc
        function l_y_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_y;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.leftWingCM(frameCurr,2) = handles.leftWingCM(frameCurr,2) + adjustVal ;
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_z_inc
        function l_z_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_z;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.leftWingCM(frameCurr,3) = handles.leftWingCM(frameCurr,3) + adjustVal ;
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_phi_inc
        function l_phi_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_phi;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.span2Hats(frameCurr,:) = ...
                    RotatePoint(handles.span2Hats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles.chord2Hats(frameCurr,:) = ...
                    RotatePoint(handles.chord2Hats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles.phi2Hats(frameCurr,:) = ...
                    RotatePoint(handles.phi2Hats(frameCurr,:),[0 0 0],[0 0 1],-adjustVal);
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_theta_inc
        function l_theta_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_theta;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.span2Hats(frameCurr,:) = ...
                    RotatePoint2(handles.span2Hats(frameCurr,:),[0 0 0],...
                    handles.phi2Hats(frameCurr,:),-1*adjustVal);
                handles.chord2Hats(frameCurr,:) = ...
                    RotatePoint2(handles.chord2Hats(frameCurr,:),[0 0 0],...
                    handles.phi2Hats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'lvec') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: l_eta_inc
        function l_eta_incButtonPushed(app, event)
            [hObject, ~, handles] = convertToGUIDECallbackArguments(app, event);
            slider = app.l_eta;
            adjustVal = 1;
            frameCurr = app.frameCurr;

            originalVal = slider.Value;
            newVal = originalVal + adjustVal;

            if newVal <= slider.Limits(2) && newVal >= slider.Limits(1)
                slider.Value = newVal;
                handles.chord2Hats(frameCurr,:) = ...
                    RotatePoint(handles.chord2Hats(frameCurr,:),[0 0 0],...
                    handles.span2Hats(frameCurr,:),-1*adjustVal);
                handles = updateChangeFlags(app, handles,'lchord') ;
                guidata(hObject, handles);
                updateDisplay(app, hObject);
            end
        end

        % Button pushed function: EZViewButton
        function EZViewButtonPushed(app, event)
            app.EZViewButton.Enable = "off";
            
            [~,~, handles] = convertToGUIDECallbackArguments(app,event);
            %disp(num2str(handles.tmsCurr))

            app.EZViewApp = correctionGUI_EZView(app, handles, app.frameCurr, app.tmsCurr);
        end

        % Close request function: figure1
        function figure1CloseRequest(app, event)
            %Close dialog box if open
            delete(app.EZViewApp)

            %Close app
            delete(app)
        end

        % Button pushed function: phi_view
        function phi_viewButtonPushed(app, event)
            [~, ~, handles] = convertToGUIDECallbackArguments(app, event);
            azview = atan2(handles.rollHats(app.frameCurr,2),handles.rollHats(app.frameCurr,1))*180/pi-90 ;
            elview = 90 ;
            view(app.main_axes,azview,elview) ;

        end

        % Button pushed function: r_theta_view
        function r_theta_viewButtonPushed(app, event)
            [~, ~, handles] = convertToGUIDECallbackArguments(app, event);
            wingBodyProduct = dot(handles.rollHats(app.frameCurr,1:2),handles.span1Hats(app.frameCurr,1:2));
            if wingBodyProduct > 0
                azview = atan2(handles.span1Hats(app.frameCurr,2),handles.span1Hats(app.frameCurr,1))*180/pi;
            else
                azview = atan2(handles.span1Hats(app.frameCurr,2),handles.span1Hats(app.frameCurr,1))*180/pi+180;
            end

            elview = 0;
            view(app.main_axes,azview,elview) ;
        end

        % Button pushed function: l_theta_view
        function l_theta_viewButtonPushed(app, event)
            [~, ~, handles] = convertToGUIDECallbackArguments(app, event);
            wingBodyProduct = dot(handles.rollHats(app.frameCurr,1:2),handles.span2Hats(app.frameCurr,1:2));
            if wingBodyProduct > 0
                azview = atan2(handles.span2Hats(app.frameCurr,2),handles.span2Hats(app.frameCurr,1))*180/pi+180;
            else
                azview = atan2(handles.span2Hats(app.frameCurr,2),handles.span2Hats(app.frameCurr,1))*180/pi;
            end

            elview = 0;
            view(app.main_axes,azview,elview);

        end

        % Button pushed function: r_eta_view
        function r_eta_viewButtonPushed(app, event)
            [~, ~, handles] = convertToGUIDECallbackArguments(app, event);
            view(app.main_axes,handles.span1Hats(app.frameCurr,:)) ;

        end

        % Button pushed function: l_eta_view
        function l_eta_viewButtonPushed(app, event)
            [~, ~, handles] = convertToGUIDECallbackArguments(app, event);
            view(app.main_axes,handles.span2Hats(app.frameCurr,:)) ;
        end

        % Button pushed function: roll_view
        function roll_viewButtonPushed(app, event)
            [~, ~, handles] = convertToGUIDECallbackArguments(app, event);
            view(app.main_axes,handles.rollHats(app.frameCurr,:)) ;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create figure1 and hide until all components are created
            app.figure1 = uifigure('Visible', 'off');
            app.figure1.Position = [680 213 1241 868];
            app.figure1.Name = 'correctionGUI_sam';
            app.figure1.CloseRequestFcn = createCallbackFcn(app, @figure1CloseRequest, true);
            app.figure1.WindowKeyPressFcn = createCallbackFcn(app, @figure1_WindowKeyPressFcn, true);
            app.figure1.HandleVisibility = 'callback';
            app.figure1.Tag = 'figure1';

            % Create uitoolbar1
            app.uitoolbar1 = uitoolbar(app.figure1);
            app.uitoolbar1.Tag = 'uitoolbar1';

            % Create ui_zoom_in
            app.ui_zoom_in = uitoggletool(app.uitoolbar1);
            app.ui_zoom_in.Tag = 'ui_zoom_in';
            app.ui_zoom_in.Tooltip = 'Zoom In';
            app.ui_zoom_in.ClickedCallback = createCallbackFcn(app, @ui_zoom_in_ClickedCallback, true);
            app.ui_zoom_in.Icon = 'ui_zoom_in_image.png';

            % Create ui_zoom_out
            app.ui_zoom_out = uitoggletool(app.uitoolbar1);
            app.ui_zoom_out.Tag = 'ui_zoom_out';
            app.ui_zoom_out.Tooltip = 'Zoom Out';
            app.ui_zoom_out.ClickedCallback = createCallbackFcn(app, @ui_zoom_out_ClickedCallback, true);
            app.ui_zoom_out.Icon = 'ui_zoom_out_image.png';

            % Create ui_rotate
            app.ui_rotate = uitoggletool(app.uitoolbar1);
            app.ui_rotate.Tag = 'ui_rotate';
            app.ui_rotate.Tooltip = 'Rotate 3D';
            app.ui_rotate.ClickedCallback = createCallbackFcn(app, @ui_rotate_ClickedCallback, true);
            app.ui_rotate.Icon = 'ui_rotate_image.png';

            % Create ui_pan
            app.ui_pan = uitoggletool(app.uitoolbar1);
            app.ui_pan.Tag = 'ui_pan';
            app.ui_pan.Tooltip = 'Pan';
            app.ui_pan.ClickedCallback = createCallbackFcn(app, @ui_pan_ClickedCallback, true);
            app.ui_pan.Icon = 'ui_pan_image.png';

            % Create ui_save
            app.ui_save = uipushtool(app.uitoolbar1);
            app.ui_save.BusyAction = 'cancel';
            app.ui_save.Interruptible = 'off';
            app.ui_save.Tag = 'ui_save';
            app.ui_save.Tooltip = 'Save Figure';
            app.ui_save.ClickedCallback = createCallbackFcn(app, @ui_save_ClickedCallback, true);
            app.ui_save.Icon = 'ui_save_image.png';

            % Create main_axes
            app.main_axes = uiaxes(app.figure1);
            app.main_axes.FontSize = 13.3333333333333;
            app.main_axes.NextPlot = 'replace';
            app.main_axes.Tag = 'main_axes';
            app.main_axes.Position = [26 76 757 729];

            % Create frame_info
            app.frame_info = uilabel(app.figure1);
            app.frame_info.Tag = 'frame_info';
            app.frame_info.BackgroundColor = [1 1 1];
            app.frame_info.HorizontalAlignment = 'center';
            app.frame_info.VerticalAlignment = 'top';
            app.frame_info.WordWrap = 'on';
            app.frame_info.FontSize = 16;
            app.frame_info.FontWeight = 'bold';
            app.frame_info.Position = [50 814 725 33];
            app.frame_info.Text = 'Frame Info';

            % Create correction_panel
            app.correction_panel = uipanel(app.figure1);
            app.correction_panel.Title = 'Correct body/wing orientation';
            app.correction_panel.Tag = 'correction_panel';
            app.correction_panel.FontSize = 10.6666666666667;
            app.correction_panel.Position = [792 15 439 643];

            % Create body_x
            app.body_x = uislider(app.correction_panel);
            app.body_x.Limits = [-80 80];
            app.body_x.MajorTicks = [];
            app.body_x.ValueChangedFcn = createCallbackFcn(app, @body_x_Callback, true);
            app.body_x.MinorTicks = [];
            app.body_x.Tag = 'body_x';
            app.body_x.Enable = 'off';
            app.body_x.FontSize = 10.6666666666667;
            app.body_x.Position = [79 597 334 3];

            % Create body_label
            app.body_label = uilabel(app.correction_panel);
            app.body_label.Tag = 'body_label';
            app.body_label.HorizontalAlignment = 'center';
            app.body_label.VerticalAlignment = 'top';
            app.body_label.WordWrap = 'on';
            app.body_label.FontSize = 13.3333333333333;
            app.body_label.FontWeight = 'bold';
            app.body_label.Position = [201 603 86 19.620895522388];
            app.body_label.Text = 'Body';

            % Create x_cm_label
            app.x_cm_label = uilabel(app.correction_panel);
            app.x_cm_label.Tag = 'x_cm_label';
            app.x_cm_label.HorizontalAlignment = 'center';
            app.x_cm_label.VerticalAlignment = 'top';
            app.x_cm_label.WordWrap = 'on';
            app.x_cm_label.FontSize = 10.6666666666667;
            app.x_cm_label.FontWeight = 'bold';
            app.x_cm_label.Position = [19 592 35.5946601941748 13.0805970149254];
            app.x_cm_label.Text = 'Xcm';

            % Create body_y
            app.body_y = uislider(app.correction_panel);
            app.body_y.Limits = [-80 80];
            app.body_y.MajorTicks = [];
            app.body_y.ValueChangedFcn = createCallbackFcn(app, @body_y_Callback, true);
            app.body_y.MinorTicks = [];
            app.body_y.Tag = 'body_y';
            app.body_y.Enable = 'off';
            app.body_y.FontSize = 10.6666666666667;
            app.body_y.Position = [79 574 334 3];

            % Create body_z
            app.body_z = uislider(app.correction_panel);
            app.body_z.Limits = [-80 80];
            app.body_z.MajorTicks = [];
            app.body_z.ValueChangedFcn = createCallbackFcn(app, @body_z_Callback, true);
            app.body_z.MinorTicks = [];
            app.body_z.Tag = 'body_z';
            app.body_z.Enable = 'off';
            app.body_z.FontSize = 10.6666666666667;
            app.body_z.Position = [79 552 334 3];

            % Create rx
            app.rx = uislider(app.correction_panel);
            app.rx.Limits = [-80 80];
            app.rx.MajorTicks = [];
            app.rx.ValueChangedFcn = createCallbackFcn(app, @rx_Callback, true);
            app.rx.MinorTicks = [];
            app.rx.Tag = 'rx';
            app.rx.Enable = 'off';
            app.rx.FontSize = 10.6666666666667;
            app.rx.Position = [79 420 334 3];

            % Create r_y
            app.r_y = uislider(app.correction_panel);
            app.r_y.Limits = [-80 80];
            app.r_y.MajorTicks = [];
            app.r_y.ValueChangedFcn = createCallbackFcn(app, @r_y_Callback, true);
            app.r_y.MinorTicks = [];
            app.r_y.Tag = 'r_y';
            app.r_y.Enable = 'off';
            app.r_y.FontSize = 10.6666666666667;
            app.r_y.Position = [79 398 334 3];

            % Create r_z
            app.r_z = uislider(app.correction_panel);
            app.r_z.Limits = [-80 80];
            app.r_z.MajorTicks = [];
            app.r_z.ValueChangedFcn = createCallbackFcn(app, @r_z_Callback, true);
            app.r_z.MinorTicks = [];
            app.r_z.Tag = 'r_z';
            app.r_z.Enable = 'off';
            app.r_z.FontSize = 10.6666666666667;
            app.r_z.Position = [79 375 334 3];

            % Create r_phi
            app.r_phi = uislider(app.correction_panel);
            app.r_phi.Limits = [-90 90];
            app.r_phi.MajorTicks = [];
            app.r_phi.ValueChangedFcn = createCallbackFcn(app, @r_phi_Callback, true);
            app.r_phi.MinorTicks = [];
            app.r_phi.Tag = 'r_phi';
            app.r_phi.Enable = 'off';
            app.r_phi.FontSize = 10.6666666666667;
            app.r_phi.Position = [79 343 334 3];

            % Create r_theta
            app.r_theta = uislider(app.correction_panel);
            app.r_theta.Limits = [-90 90];
            app.r_theta.MajorTicks = [];
            app.r_theta.ValueChangedFcn = createCallbackFcn(app, @r_theta_Callback, true);
            app.r_theta.MinorTicks = [];
            app.r_theta.Tag = 'r_theta';
            app.r_theta.Enable = 'off';
            app.r_theta.FontSize = 10.6666666666667;
            app.r_theta.Position = [79 317 334 3];

            % Create r_eta
            app.r_eta = uislider(app.correction_panel);
            app.r_eta.Limits = [-180 180];
            app.r_eta.MajorTicks = [];
            app.r_eta.ValueChangedFcn = createCallbackFcn(app, @r_eta_Callback, true);
            app.r_eta.MinorTicks = [];
            app.r_eta.Tag = 'r_eta';
            app.r_eta.Enable = 'off';
            app.r_eta.FontSize = 10.6666666666667;
            app.r_eta.Position = [79 291 334 3];

            % Create l_x
            app.l_x = uislider(app.correction_panel);
            app.l_x.Limits = [-80 80];
            app.l_x.MajorTicks = [];
            app.l_x.ValueChangedFcn = createCallbackFcn(app, @l_x_Callback, true);
            app.l_x.MinorTicks = [];
            app.l_x.Tag = 'l_x';
            app.l_x.Enable = 'off';
            app.l_x.FontSize = 10.6666666666667;
            app.l_x.Position = [79 234 334 3];

            % Create l_y
            app.l_y = uislider(app.correction_panel);
            app.l_y.Limits = [-80 80];
            app.l_y.MajorTicks = [];
            app.l_y.ValueChangedFcn = createCallbackFcn(app, @l_y_Callback, true);
            app.l_y.MinorTicks = [];
            app.l_y.Tag = 'l_y';
            app.l_y.Enable = 'off';
            app.l_y.FontSize = 10.6666666666667;
            app.l_y.Position = [79 212 334 3];

            % Create l_z
            app.l_z = uislider(app.correction_panel);
            app.l_z.Limits = [-80 80];
            app.l_z.MajorTicks = [];
            app.l_z.ValueChangedFcn = createCallbackFcn(app, @l_z_Callback, true);
            app.l_z.MinorTicks = [];
            app.l_z.Tag = 'l_z';
            app.l_z.Enable = 'off';
            app.l_z.FontSize = 10.6666666666667;
            app.l_z.Position = [79 189 334 3];

            % Create l_phi
            app.l_phi = uislider(app.correction_panel);
            app.l_phi.Limits = [-90 90];
            app.l_phi.MajorTicks = [];
            app.l_phi.ValueChangedFcn = createCallbackFcn(app, @l_phi_Callback, true);
            app.l_phi.MinorTicks = [];
            app.l_phi.Tag = 'l_phi';
            app.l_phi.Enable = 'off';
            app.l_phi.FontSize = 10.6666666666667;
            app.l_phi.Position = [79 158 334 3];

            % Create l_theta
            app.l_theta = uislider(app.correction_panel);
            app.l_theta.Limits = [-90 90];
            app.l_theta.MajorTicks = [];
            app.l_theta.ValueChangedFcn = createCallbackFcn(app, @l_theta_Callback, true);
            app.l_theta.MinorTicks = [];
            app.l_theta.Tag = 'l_theta';
            app.l_theta.Enable = 'off';
            app.l_theta.FontSize = 10.6666666666667;
            app.l_theta.Position = [79 132 334 3];

            % Create l_eta
            app.l_eta = uislider(app.correction_panel);
            app.l_eta.Limits = [-180 180];
            app.l_eta.MajorTicks = [];
            app.l_eta.ValueChangedFcn = createCallbackFcn(app, @l_eta_Callback, true);
            app.l_eta.MinorTicks = [];
            app.l_eta.Tag = 'l_eta';
            app.l_eta.Enable = 'off';
            app.l_eta.FontSize = 10.6666666666667;
            app.l_eta.Position = [79 107 334 3];

            % Create right_wing_label
            app.right_wing_label = uilabel(app.correction_panel);
            app.right_wing_label.Tag = 'right_wing_label';
            app.right_wing_label.HorizontalAlignment = 'center';
            app.right_wing_label.VerticalAlignment = 'top';
            app.right_wing_label.WordWrap = 'on';
            app.right_wing_label.FontSize = 13.3333333333333;
            app.right_wing_label.FontWeight = 'bold';
            app.right_wing_label.Position = [199 436 90 16.8179104477612];
            app.right_wing_label.Text = 'Right Wing';

            % Create left_wing_label
            app.left_wing_label = uilabel(app.correction_panel);
            app.left_wing_label.Tag = 'left_wing_label';
            app.left_wing_label.HorizontalAlignment = 'center';
            app.left_wing_label.VerticalAlignment = 'top';
            app.left_wing_label.WordWrap = 'on';
            app.left_wing_label.FontSize = 13.3333333333333;
            app.left_wing_label.FontWeight = 'bold';
            app.left_wing_label.Position = [207 250 71 19.6208955223881];
            app.left_wing_label.Text = 'Left Wing';

            % Create body_phi
            app.body_phi = uislider(app.correction_panel);
            app.body_phi.Limits = [-180 180];
            app.body_phi.MajorTicks = [];
            app.body_phi.ValueChangedFcn = createCallbackFcn(app, @body_phi_Callback, true);
            app.body_phi.MinorTicks = [];
            app.body_phi.Tag = 'body_phi';
            app.body_phi.Enable = 'off';
            app.body_phi.FontSize = 10.6666666666667;
            app.body_phi.Position = [79 521 334 3];

            % Create body_pitch
            app.body_pitch = uislider(app.correction_panel);
            app.body_pitch.Limits = [-180 180];
            app.body_pitch.MajorTicks = [];
            app.body_pitch.ValueChangedFcn = createCallbackFcn(app, @body_pitch_Callback, true);
            app.body_pitch.MinorTicks = [];
            app.body_pitch.Tag = 'body_pitch';
            app.body_pitch.Enable = 'off';
            app.body_pitch.FontSize = 10.6666666666667;
            app.body_pitch.Position = [79 498 334 3];

            % Create body_roll
            app.body_roll = uislider(app.correction_panel);
            app.body_roll.Limits = [-180 180];
            app.body_roll.MajorTicks = [];
            app.body_roll.ValueChangedFcn = createCallbackFcn(app, @body_roll_Callback, true);
            app.body_roll.MinorTicks = [];
            app.body_roll.Tag = 'body_roll';
            app.body_roll.Enable = 'off';
            app.body_roll.FontSize = 10.6666666666667;
            app.body_roll.Position = [79 476 334 3];

            % Create ignoreFrame
            app.ignoreFrame = uibutton(app.correction_panel, 'state');
            app.ignoreFrame.ValueChangedFcn = createCallbackFcn(app, @ignoreFrame_Callback, true);
            app.ignoreFrame.Tag = 'ignoreFrame';
            app.ignoreFrame.Enable = 'off';
            app.ignoreFrame.Text = 'Ignore Frame';
            app.ignoreFrame.FontSize = 10.6666666666667;
            app.ignoreFrame.Position = [20 47 126 26.8285714285714];

            % Create swap_wings
            app.swap_wings = uibutton(app.correction_panel, 'push');
            app.swap_wings.ButtonPushedFcn = createCallbackFcn(app, @swap_wings_Callback, true);
            app.swap_wings.Tag = 'swap_wings';
            app.swap_wings.FontSize = 10.6666666666667;
            app.swap_wings.Enable = 'off';
            app.swap_wings.Position = [292 47 126 26.8285714285714];
            app.swap_wings.Text = 'Swap Wings';

            % Create y_cm_label
            app.y_cm_label = uilabel(app.correction_panel);
            app.y_cm_label.Tag = 'y_cm_label';
            app.y_cm_label.HorizontalAlignment = 'center';
            app.y_cm_label.VerticalAlignment = 'top';
            app.y_cm_label.WordWrap = 'on';
            app.y_cm_label.FontSize = 10.6666666666667;
            app.y_cm_label.FontWeight = 'bold';
            app.y_cm_label.Position = [19 569 35.5946601941748 13.0805970149254];
            app.y_cm_label.Text = 'Ycm';

            % Create z_cm_label
            app.z_cm_label = uilabel(app.correction_panel);
            app.z_cm_label.Tag = 'z_cm_label';
            app.z_cm_label.HorizontalAlignment = 'center';
            app.z_cm_label.VerticalAlignment = 'top';
            app.z_cm_label.WordWrap = 'on';
            app.z_cm_label.FontSize = 10.6666666666667;
            app.z_cm_label.FontWeight = 'bold';
            app.z_cm_label.Position = [19 547 35.5946601941748 13.0805970149254];
            app.z_cm_label.Text = 'Zcm';

            % Create yaw_label
            app.yaw_label = uilabel(app.correction_panel);
            app.yaw_label.Tag = 'yaw_label';
            app.yaw_label.HorizontalAlignment = 'center';
            app.yaw_label.VerticalAlignment = 'top';
            app.yaw_label.WordWrap = 'on';
            app.yaw_label.FontSize = 10.6666666666667;
            app.yaw_label.FontWeight = 'bold';
            app.yaw_label.Position = [19 516 35.5946601941748 13.0805970149254];
            app.yaw_label.Text = 'Yaw';

            % Create pitch_label
            app.pitch_label = uilabel(app.correction_panel);
            app.pitch_label.Tag = 'pitch_label';
            app.pitch_label.HorizontalAlignment = 'center';
            app.pitch_label.VerticalAlignment = 'top';
            app.pitch_label.WordWrap = 'on';
            app.pitch_label.FontSize = 10.6666666666667;
            app.pitch_label.FontWeight = 'bold';
            app.pitch_label.Position = [19 493 35.5946601941748 13.0805970149254];
            app.pitch_label.Text = 'Pitch';

            % Create roll_label
            app.roll_label = uilabel(app.correction_panel);
            app.roll_label.Tag = 'roll_label';
            app.roll_label.HorizontalAlignment = 'center';
            app.roll_label.VerticalAlignment = 'top';
            app.roll_label.WordWrap = 'on';
            app.roll_label.FontSize = 10.6666666666667;
            app.roll_label.FontWeight = 'bold';
            app.roll_label.Position = [19 471 35.5946601941748 13.0805970149254];
            app.roll_label.Text = 'Roll';

            % Create r_x_label
            app.r_x_label = uilabel(app.correction_panel);
            app.r_x_label.Tag = 'r_x_label';
            app.r_x_label.HorizontalAlignment = 'center';
            app.r_x_label.VerticalAlignment = 'top';
            app.r_x_label.WordWrap = 'on';
            app.r_x_label.FontSize = 10.6666666666667;
            app.r_x_label.FontWeight = 'bold';
            app.r_x_label.FontColor = [0.64 0.08 0.18];
            app.r_x_label.Position = [21 415 36 13.0805970149254];
            app.r_x_label.Text = 'R X';

            % Create r_y_label
            app.r_y_label = uilabel(app.correction_panel);
            app.r_y_label.Tag = 'r_y_label';
            app.r_y_label.HorizontalAlignment = 'center';
            app.r_y_label.VerticalAlignment = 'top';
            app.r_y_label.WordWrap = 'on';
            app.r_y_label.FontSize = 10.6666666666667;
            app.r_y_label.FontWeight = 'bold';
            app.r_y_label.FontColor = [0.64 0.08 0.18];
            app.r_y_label.Position = [21 393 36 13.0805970149254];
            app.r_y_label.Text = 'R Y';

            % Create r_z_label
            app.r_z_label = uilabel(app.correction_panel);
            app.r_z_label.Tag = 'r_z_label';
            app.r_z_label.HorizontalAlignment = 'center';
            app.r_z_label.VerticalAlignment = 'top';
            app.r_z_label.WordWrap = 'on';
            app.r_z_label.FontSize = 10.6666666666667;
            app.r_z_label.FontWeight = 'bold';
            app.r_z_label.FontColor = [0.64 0.08 0.18];
            app.r_z_label.Position = [21 370 36 13.0805970149254];
            app.r_z_label.Text = 'R Z';

            % Create l_x_label
            app.l_x_label = uilabel(app.correction_panel);
            app.l_x_label.Tag = 'l_x_label';
            app.l_x_label.HorizontalAlignment = 'center';
            app.l_x_label.VerticalAlignment = 'top';
            app.l_x_label.WordWrap = 'on';
            app.l_x_label.FontSize = 10.6666666666667;
            app.l_x_label.FontWeight = 'bold';
            app.l_x_label.FontColor = [0 0.45 0.74];
            app.l_x_label.Position = [22 229 36 13.0805970149254];
            app.l_x_label.Text = 'L X';

            % Create l_y_label
            app.l_y_label = uilabel(app.correction_panel);
            app.l_y_label.Tag = 'l_y_label';
            app.l_y_label.HorizontalAlignment = 'center';
            app.l_y_label.VerticalAlignment = 'top';
            app.l_y_label.WordWrap = 'on';
            app.l_y_label.FontSize = 10.6666666666667;
            app.l_y_label.FontWeight = 'bold';
            app.l_y_label.FontColor = [0 0.45 0.74];
            app.l_y_label.Position = [22 207 36 13.0805970149254];
            app.l_y_label.Text = 'L Y';

            % Create l_z_label
            app.l_z_label = uilabel(app.correction_panel);
            app.l_z_label.Tag = 'l_z_label';
            app.l_z_label.HorizontalAlignment = 'center';
            app.l_z_label.VerticalAlignment = 'top';
            app.l_z_label.WordWrap = 'on';
            app.l_z_label.FontSize = 10.6666666666667;
            app.l_z_label.FontWeight = 'bold';
            app.l_z_label.FontColor = [0 0.45 0.74];
            app.l_z_label.Position = [22 184 36 13.0805970149254];
            app.l_z_label.Text = 'L Z';

            % Create cluster_wings
            app.cluster_wings = uibutton(app.correction_panel, 'push');
            app.cluster_wings.ButtonPushedFcn = createCallbackFcn(app, @cluster_wings_Callback, true);
            app.cluster_wings.Tag = 'cluster_wings';
            app.cluster_wings.FontSize = 10.6666666666667;
            app.cluster_wings.Enable = 'off';
            app.cluster_wings.Position = [156 47 126 26.8285714285714];
            app.cluster_wings.Text = 'Cluster Wings';

            % Create save_manualCorrRangeMS
            app.save_manualCorrRangeMS = uibutton(app.correction_panel, 'push');
            app.save_manualCorrRangeMS.ButtonPushedFcn = createCallbackFcn(app, @save_manualCorrRangeMS_Callback, true);
            app.save_manualCorrRangeMS.Tag = 'save_manualCorrRangeMS';
            app.save_manualCorrRangeMS.FontSize = 10.6666666666667;
            app.save_manualCorrRangeMS.Enable = 'off';
            app.save_manualCorrRangeMS.Position = [292 11 126 27];
            app.save_manualCorrRangeMS.Text = 'Enter Manual Corr Range';

            % Create save_roll
            app.save_roll = uibutton(app.correction_panel, 'state');
            app.save_roll.ValueChangedFcn = createCallbackFcn(app, @save_roll_Callback, true);
            app.save_roll.Tag = 'save_roll';
            app.save_roll.Enable = 'off';
            app.save_roll.Text = 'Save Roll';
            app.save_roll.FontSize = 10.6666666666667;
            app.save_roll.Position = [22 11 124 27];

            % Create body_x_dec
            app.body_x_dec = uibutton(app.correction_panel, 'push');
            app.body_x_dec.ButtonPushedFcn = createCallbackFcn(app, @body_x_decButtonPushed, true);
            app.body_x_dec.HorizontalAlignment = 'left';
            app.body_x_dec.Enable = 'off';
            app.body_x_dec.Position = [62 587 18 23];
            app.body_x_dec.Text = '◀';

            % Create body_x_inc
            app.body_x_inc = uibutton(app.correction_panel, 'push');
            app.body_x_inc.ButtonPushedFcn = createCallbackFcn(app, @body_x_incButtonPushed, true);
            app.body_x_inc.Enable = 'off';
            app.body_x_inc.Position = [412 587 18 23];
            app.body_x_inc.Text = '▶';

            % Create body_y_dec
            app.body_y_dec = uibutton(app.correction_panel, 'push');
            app.body_y_dec.ButtonPushedFcn = createCallbackFcn(app, @body_y_decButtonPushed, true);
            app.body_y_dec.HorizontalAlignment = 'left';
            app.body_y_dec.Enable = 'off';
            app.body_y_dec.Position = [62 564 18 23];
            app.body_y_dec.Text = '◀';

            % Create body_y_inc
            app.body_y_inc = uibutton(app.correction_panel, 'push');
            app.body_y_inc.ButtonPushedFcn = createCallbackFcn(app, @body_y_incButtonPushed, true);
            app.body_y_inc.Enable = 'off';
            app.body_y_inc.Position = [412 564 18 23];
            app.body_y_inc.Text = '▶';

            % Create body_z_dec
            app.body_z_dec = uibutton(app.correction_panel, 'push');
            app.body_z_dec.ButtonPushedFcn = createCallbackFcn(app, @body_z_decButtonPushed, true);
            app.body_z_dec.HorizontalAlignment = 'left';
            app.body_z_dec.Enable = 'off';
            app.body_z_dec.Position = [62 542 18 23];
            app.body_z_dec.Text = '◀';

            % Create body_phi_dec
            app.body_phi_dec = uibutton(app.correction_panel, 'push');
            app.body_phi_dec.ButtonPushedFcn = createCallbackFcn(app, @body_phi_decButtonPushed, true);
            app.body_phi_dec.HorizontalAlignment = 'left';
            app.body_phi_dec.Enable = 'off';
            app.body_phi_dec.Position = [62 511 18 23];
            app.body_phi_dec.Text = '◀';

            % Create body_pitch_dec
            app.body_pitch_dec = uibutton(app.correction_panel, 'push');
            app.body_pitch_dec.ButtonPushedFcn = createCallbackFcn(app, @body_pitch_decButtonPushed, true);
            app.body_pitch_dec.HorizontalAlignment = 'left';
            app.body_pitch_dec.Enable = 'off';
            app.body_pitch_dec.Position = [62 488 18 23];
            app.body_pitch_dec.Text = '◀';

            % Create body_roll_dec
            app.body_roll_dec = uibutton(app.correction_panel, 'push');
            app.body_roll_dec.ButtonPushedFcn = createCallbackFcn(app, @body_roll_decButtonPushed, true);
            app.body_roll_dec.HorizontalAlignment = 'left';
            app.body_roll_dec.Enable = 'off';
            app.body_roll_dec.Position = [62 466 18 23];
            app.body_roll_dec.Text = '◀';

            % Create rx_dec
            app.rx_dec = uibutton(app.correction_panel, 'push');
            app.rx_dec.ButtonPushedFcn = createCallbackFcn(app, @rx_decButtonPushed, true);
            app.rx_dec.HorizontalAlignment = 'left';
            app.rx_dec.Enable = 'off';
            app.rx_dec.Position = [62 410 18 23];
            app.rx_dec.Text = '◀';

            % Create r_y_dec
            app.r_y_dec = uibutton(app.correction_panel, 'push');
            app.r_y_dec.ButtonPushedFcn = createCallbackFcn(app, @r_y_decButtonPushed, true);
            app.r_y_dec.HorizontalAlignment = 'left';
            app.r_y_dec.Enable = 'off';
            app.r_y_dec.Position = [62 388 18 23];
            app.r_y_dec.Text = '◀';

            % Create r_z_dec
            app.r_z_dec = uibutton(app.correction_panel, 'push');
            app.r_z_dec.ButtonPushedFcn = createCallbackFcn(app, @r_z_decButtonPushed, true);
            app.r_z_dec.HorizontalAlignment = 'left';
            app.r_z_dec.Enable = 'off';
            app.r_z_dec.Position = [62 365 18 23];
            app.r_z_dec.Text = '◀';

            % Create r_phi_dec
            app.r_phi_dec = uibutton(app.correction_panel, 'push');
            app.r_phi_dec.ButtonPushedFcn = createCallbackFcn(app, @r_phi_decButtonPushed, true);
            app.r_phi_dec.HorizontalAlignment = 'left';
            app.r_phi_dec.Enable = 'off';
            app.r_phi_dec.Position = [62 333 18 23];
            app.r_phi_dec.Text = '◀';

            % Create r_theta_dec
            app.r_theta_dec = uibutton(app.correction_panel, 'push');
            app.r_theta_dec.ButtonPushedFcn = createCallbackFcn(app, @r_theta_decButtonPushed, true);
            app.r_theta_dec.HorizontalAlignment = 'left';
            app.r_theta_dec.Enable = 'off';
            app.r_theta_dec.Position = [62 307 18 23];
            app.r_theta_dec.Text = '◀';

            % Create r_eta_dec
            app.r_eta_dec = uibutton(app.correction_panel, 'push');
            app.r_eta_dec.ButtonPushedFcn = createCallbackFcn(app, @r_eta_decButtonPushed, true);
            app.r_eta_dec.HorizontalAlignment = 'left';
            app.r_eta_dec.Enable = 'off';
            app.r_eta_dec.Position = [62 281 18 23];
            app.r_eta_dec.Text = '◀';

            % Create l_x_dec
            app.l_x_dec = uibutton(app.correction_panel, 'push');
            app.l_x_dec.ButtonPushedFcn = createCallbackFcn(app, @l_x_decButtonPushed, true);
            app.l_x_dec.HorizontalAlignment = 'left';
            app.l_x_dec.Enable = 'off';
            app.l_x_dec.Position = [62 224 18 23];
            app.l_x_dec.Text = '◀';

            % Create l_y_dec
            app.l_y_dec = uibutton(app.correction_panel, 'push');
            app.l_y_dec.ButtonPushedFcn = createCallbackFcn(app, @l_y_decButtonPushed, true);
            app.l_y_dec.HorizontalAlignment = 'left';
            app.l_y_dec.Enable = 'off';
            app.l_y_dec.Position = [62 202 18 23];
            app.l_y_dec.Text = '◀';

            % Create l_z_dec
            app.l_z_dec = uibutton(app.correction_panel, 'push');
            app.l_z_dec.ButtonPushedFcn = createCallbackFcn(app, @l_z_decButtonPushed, true);
            app.l_z_dec.HorizontalAlignment = 'left';
            app.l_z_dec.Enable = 'off';
            app.l_z_dec.Position = [62 179 18 23];
            app.l_z_dec.Text = '◀';

            % Create l_phi_dec
            app.l_phi_dec = uibutton(app.correction_panel, 'push');
            app.l_phi_dec.ButtonPushedFcn = createCallbackFcn(app, @l_phi_decButtonPushed, true);
            app.l_phi_dec.HorizontalAlignment = 'left';
            app.l_phi_dec.Enable = 'off';
            app.l_phi_dec.Position = [62 148 18 23];
            app.l_phi_dec.Text = '◀';

            % Create l_theta_dec
            app.l_theta_dec = uibutton(app.correction_panel, 'push');
            app.l_theta_dec.ButtonPushedFcn = createCallbackFcn(app, @l_theta_decButtonPushed, true);
            app.l_theta_dec.HorizontalAlignment = 'left';
            app.l_theta_dec.Enable = 'off';
            app.l_theta_dec.Position = [62 122 18 23];
            app.l_theta_dec.Text = '◀';

            % Create l_eta_dec
            app.l_eta_dec = uibutton(app.correction_panel, 'push');
            app.l_eta_dec.ButtonPushedFcn = createCallbackFcn(app, @l_eta_decButtonPushed, true);
            app.l_eta_dec.HorizontalAlignment = 'left';
            app.l_eta_dec.Enable = 'off';
            app.l_eta_dec.Position = [62 97 18 23];
            app.l_eta_dec.Text = '◀';

            % Create body_z_inc
            app.body_z_inc = uibutton(app.correction_panel, 'push');
            app.body_z_inc.ButtonPushedFcn = createCallbackFcn(app, @body_z_incButtonPushed, true);
            app.body_z_inc.Enable = 'off';
            app.body_z_inc.Position = [412 542 18 23];
            app.body_z_inc.Text = '▶';

            % Create body_phi_inc
            app.body_phi_inc = uibutton(app.correction_panel, 'push');
            app.body_phi_inc.ButtonPushedFcn = createCallbackFcn(app, @body_phi_incButtonPushed, true);
            app.body_phi_inc.Enable = 'off';
            app.body_phi_inc.Position = [412 511 18 23];
            app.body_phi_inc.Text = '▶';

            % Create body_pitch_inc
            app.body_pitch_inc = uibutton(app.correction_panel, 'push');
            app.body_pitch_inc.ButtonPushedFcn = createCallbackFcn(app, @body_pitch_incButtonPushed, true);
            app.body_pitch_inc.Enable = 'off';
            app.body_pitch_inc.Position = [412 488 18 23];
            app.body_pitch_inc.Text = '▶';

            % Create body_roll_inc
            app.body_roll_inc = uibutton(app.correction_panel, 'push');
            app.body_roll_inc.ButtonPushedFcn = createCallbackFcn(app, @body_roll_incButtonPushed, true);
            app.body_roll_inc.Enable = 'off';
            app.body_roll_inc.Position = [412 466 18 23];
            app.body_roll_inc.Text = '▶';

            % Create rx_inc
            app.rx_inc = uibutton(app.correction_panel, 'push');
            app.rx_inc.ButtonPushedFcn = createCallbackFcn(app, @rx_incButtonPushed, true);
            app.rx_inc.Enable = 'off';
            app.rx_inc.Position = [412 410 18 23];
            app.rx_inc.Text = '▶';

            % Create r_y_inc
            app.r_y_inc = uibutton(app.correction_panel, 'push');
            app.r_y_inc.ButtonPushedFcn = createCallbackFcn(app, @r_y_incButtonPushed, true);
            app.r_y_inc.Enable = 'off';
            app.r_y_inc.Position = [412 388 18 23];
            app.r_y_inc.Text = '▶';

            % Create r_z_inc
            app.r_z_inc = uibutton(app.correction_panel, 'push');
            app.r_z_inc.ButtonPushedFcn = createCallbackFcn(app, @r_z_incButtonPushed, true);
            app.r_z_inc.Enable = 'off';
            app.r_z_inc.Position = [412 365 18 23];
            app.r_z_inc.Text = '▶';

            % Create r_phi_inc
            app.r_phi_inc = uibutton(app.correction_panel, 'push');
            app.r_phi_inc.ButtonPushedFcn = createCallbackFcn(app, @r_phi_incButtonPushed, true);
            app.r_phi_inc.Enable = 'off';
            app.r_phi_inc.Position = [412 333 18 23];
            app.r_phi_inc.Text = '▶';

            % Create r_theta_inc
            app.r_theta_inc = uibutton(app.correction_panel, 'push');
            app.r_theta_inc.ButtonPushedFcn = createCallbackFcn(app, @r_theta_incButtonPushed, true);
            app.r_theta_inc.Enable = 'off';
            app.r_theta_inc.Position = [412 307 18 23];
            app.r_theta_inc.Text = '▶';

            % Create r_eta_inc
            app.r_eta_inc = uibutton(app.correction_panel, 'push');
            app.r_eta_inc.ButtonPushedFcn = createCallbackFcn(app, @r_eta_incButtonPushed, true);
            app.r_eta_inc.Enable = 'off';
            app.r_eta_inc.Position = [412 281 18 23];
            app.r_eta_inc.Text = '▶';

            % Create l_x_inc
            app.l_x_inc = uibutton(app.correction_panel, 'push');
            app.l_x_inc.ButtonPushedFcn = createCallbackFcn(app, @l_x_incButtonPushed, true);
            app.l_x_inc.Enable = 'off';
            app.l_x_inc.Position = [412 224 18 23];
            app.l_x_inc.Text = '▶';

            % Create l_y_inc
            app.l_y_inc = uibutton(app.correction_panel, 'push');
            app.l_y_inc.ButtonPushedFcn = createCallbackFcn(app, @l_y_incButtonPushed, true);
            app.l_y_inc.Enable = 'off';
            app.l_y_inc.Position = [412 202 18 23];
            app.l_y_inc.Text = '▶';

            % Create l_z_inc
            app.l_z_inc = uibutton(app.correction_panel, 'push');
            app.l_z_inc.ButtonPushedFcn = createCallbackFcn(app, @l_z_incButtonPushed, true);
            app.l_z_inc.Enable = 'off';
            app.l_z_inc.Position = [412 179 18 23];
            app.l_z_inc.Text = '▶';

            % Create l_phi_inc
            app.l_phi_inc = uibutton(app.correction_panel, 'push');
            app.l_phi_inc.ButtonPushedFcn = createCallbackFcn(app, @l_phi_incButtonPushed, true);
            app.l_phi_inc.Enable = 'off';
            app.l_phi_inc.Position = [412 148 18 23];
            app.l_phi_inc.Text = '▶';

            % Create l_theta_inc
            app.l_theta_inc = uibutton(app.correction_panel, 'push');
            app.l_theta_inc.ButtonPushedFcn = createCallbackFcn(app, @l_theta_incButtonPushed, true);
            app.l_theta_inc.Enable = 'off';
            app.l_theta_inc.Position = [412 122 18 23];
            app.l_theta_inc.Text = '▶';

            % Create l_eta_inc
            app.l_eta_inc = uibutton(app.correction_panel, 'push');
            app.l_eta_inc.ButtonPushedFcn = createCallbackFcn(app, @l_eta_incButtonPushed, true);
            app.l_eta_inc.Enable = 'off';
            app.l_eta_inc.Position = [412 96 18 23];
            app.l_eta_inc.Text = '▶';

            % Create EZViewButton
            app.EZViewButton = uibutton(app.correction_panel, 'push');
            app.EZViewButton.ButtonPushedFcn = createCallbackFcn(app, @EZViewButtonPushed, true);
            app.EZViewButton.Enable = 'off';
            app.EZViewButton.Position = [156 11 126 27];
            app.EZViewButton.Text = 'EZ View';

            % Create r_theta_view
            app.r_theta_view = uibutton(app.correction_panel, 'push');
            app.r_theta_view.ButtonPushedFcn = createCallbackFcn(app, @r_theta_viewButtonPushed, true);
            app.r_theta_view.FontColor = [0.6392 0.0784 0.1804];
            app.r_theta_view.Enable = 'off';
            app.r_theta_view.Position = [5 305 57 27];
            app.r_theta_view.Text = 'R Theta';

            % Create l_eta_view
            app.l_eta_view = uibutton(app.correction_panel, 'push');
            app.l_eta_view.ButtonPushedFcn = createCallbackFcn(app, @l_eta_viewButtonPushed, true);
            app.l_eta_view.FontColor = [0.6392 0.0784 0.1804];
            app.l_eta_view.Enable = 'off';
            app.l_eta_view.Position = [5 279 57 27];
            app.l_eta_view.Text = 'L Eta';

            % Create phi_view_2
            app.phi_view_2 = uibutton(app.correction_panel, 'push');
            app.phi_view_2.FontColor = [0.6392 0.0784 0.1804];
            app.phi_view_2.Enable = 'off';
            app.phi_view_2.Position = [5 331 57 27];
            app.phi_view_2.Text = 'R Phi';

            % Create l_theta_view
            app.l_theta_view = uibutton(app.correction_panel, 'push');
            app.l_theta_view.ButtonPushedFcn = createCallbackFcn(app, @l_theta_viewButtonPushed, true);
            app.l_theta_view.FontColor = [0 0.451 0.7412];
            app.l_theta_view.Enable = 'off';
            app.l_theta_view.Position = [5 120 57 27];
            app.l_theta_view.Text = 'L Theta';

            % Create phi_view
            app.phi_view = uibutton(app.correction_panel, 'push');
            app.phi_view.ButtonPushedFcn = createCallbackFcn(app, @phi_viewButtonPushed, true);
            app.phi_view.FontColor = [0 0.451 0.7412];
            app.phi_view.Enable = 'off';
            app.phi_view.Position = [5 146 57 27];
            app.phi_view.Text = 'L Phi';

            % Create r_eta_view
            app.r_eta_view = uibutton(app.correction_panel, 'push');
            app.r_eta_view.ButtonPushedFcn = createCallbackFcn(app, @r_eta_viewButtonPushed, true);
            app.r_eta_view.FontColor = [0 0.451 0.7412];
            app.r_eta_view.Enable = 'off';
            app.r_eta_view.Position = [5 94 57 27];
            app.r_eta_view.Text = 'R Eta';

            % Create data_dir
            app.data_dir = uilistbox(app.figure1);
            app.data_dir.Items = {'(import data files into this list box)'};
            app.data_dir.ValueChangedFcn = createCallbackFcn(app, @data_dir_Callback, true);
            app.data_dir.Tag = 'data_dir';
            app.data_dir.FontSize = 10.6666666666667;
            app.data_dir.Position = [807 699 416 124];
            app.data_dir.Value = '(import data files into this list box)';

            % Create open_data_dir
            app.open_data_dir = uibutton(app.figure1, 'push');
            app.open_data_dir.ButtonPushedFcn = createCallbackFcn(app, @open_data_dir_Callback, true);
            app.open_data_dir.Tag = 'open_data_dir';
            app.open_data_dir.FontSize = 10.6666666666667;
            app.open_data_dir.Position = [808 827 415 30];
            app.open_data_dir.Text = 'Open Data Directory';

            % Create load_data
            app.load_data = uibutton(app.figure1, 'push');
            app.load_data.ButtonPushedFcn = createCallbackFcn(app, @load_data_Callback, true);
            app.load_data.Tag = 'load_data';
            app.load_data.FontSize = 10.6666666666667;
            app.load_data.Position = [808 663 202 32];
            app.load_data.Text = 'Load Data';

            % Create bback
            app.bback = uibutton(app.figure1, 'push');
            app.bback.ButtonPushedFcn = createCallbackFcn(app, @bback_Callback, true);
            app.bback.Tag = 'bback';
            app.bback.FontSize = 10.6666666666667;
            app.bback.Enable = 'off';
            app.bback.Position = [66 24 101 29];
            app.bback.Text = '<<';

            % Create back
            app.back = uibutton(app.figure1, 'push');
            app.back.ButtonPushedFcn = createCallbackFcn(app, @back_Callback, true);
            app.back.Tag = 'back';
            app.back.FontSize = 10.6666666666667;
            app.back.Enable = 'off';
            app.back.Position = [194 24 101 29];
            app.back.Text = '<';

            % Create fwd
            app.fwd = uibutton(app.figure1, 'push');
            app.fwd.ButtonPushedFcn = createCallbackFcn(app, @fwd_Callback, true);
            app.fwd.Tag = 'fwd';
            app.fwd.FontSize = 10.6666666666667;
            app.fwd.Enable = 'off';
            app.fwd.Position = [521 24 101 29];
            app.fwd.Text = '>';

            % Create ffwd
            app.ffwd = uibutton(app.figure1, 'push');
            app.ffwd.ButtonPushedFcn = createCallbackFcn(app, @ffwd_Callback, true);
            app.ffwd.Tag = 'ffwd';
            app.ffwd.FontSize = 10.6666666666667;
            app.ffwd.Enable = 'off';
            app.ffwd.Position = [644 24 101 29];
            app.ffwd.Text = '>>';

            % Create clear_data
            app.clear_data = uibutton(app.figure1, 'push');
            app.clear_data.ButtonPushedFcn = createCallbackFcn(app, @clear_data_Callback, true);
            app.clear_data.Tag = 'clear_data';
            app.clear_data.FontSize = 10.6666666666667;
            app.clear_data.Position = [1020 663 202 32];
            app.clear_data.Text = 'Clear Data';

            % Create roll_view
            app.roll_view = uibutton(app.figure1, 'push');
            app.roll_view.ButtonPushedFcn = createCallbackFcn(app, @roll_viewButtonPushed, true);
            app.roll_view.Enable = 'off';
            app.roll_view.Position = [387 25 57 27];
            app.roll_view.Text = 'Roll';

            % Show the figure after all components are created
            app.figure1.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = correctionGUI_sam_App_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.figure1)

                % Execute the startup function
                runStartupFcn(app, @(app)correctionGUI_sam_OpeningFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.figure1)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.figure1)
        end
    end
end