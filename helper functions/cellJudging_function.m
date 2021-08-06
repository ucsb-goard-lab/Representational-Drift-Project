function cellJudging_function(maps, currCell, cellInfo, fig, stopFlag)
% controls interactive gui for visual inspection of ROIs. Handled in an automated fashion by the preprocessing pipeline
    % initialize on first iteration
    if nargin < 4
        stopFlag = 0;
        fig = figure;
        scrn_size = get(0,'ScreenSize');
        set(fig,'Position',[50 50 scrn_size(3)-100 scrn_size(4)-150]);
        set(fig,'color',[1 1 1]);
    end
    if nargin < 2
        currCell = 1;
        cellInfo = struct();
    end
    if nargin < 1
        maps_filename = uigetfile('.mat');
        maps = importdata(maps_filename);
    end

    % stop function if button is pressed
    if stopFlag==1
        close
        clear global data
        return
    end

    numCells = size(maps.NatMov.centers, 1);
    num_sessions = size(maps.NatMov.avgproj_mat, 3);

    if isempty(fieldnames(cellInfo))            %if cellInfo not provided, initialize it
        cellInfo.quality = zeros(1, numCells);
        cellInfo.presence = zeros(numCells, num_sessions);
    end
    cellInfo.presence(currCell, :) = 1;     %default to 1 (present)
    cellInfo.lastCell = currCell;

    padding = maps.NatMov.padding;
    yPixels = maps.NatMov.yPixels;
    xPixels = maps.NatMov.xPixels;
    centers = maps.NatMov.centers-padding;
    curr_full_avgproj = maps.NatMov.avgproj_mat(padding+1:padding+yPixels, padding+1:padding+xPixels, 1);     %average projection from natmov on day 0
    curr_natmov_avgprojs = maps.NatMov.cell_avgproj_mat(:, :, :, currCell);
    curr_natmov_actmaps = maps.NatMov.cell_actmap_mat(:, :, :, currCell);
    curr_pdg_avgprojs = maps.PDG.cell_avgproj_mat(:, :, :, currCell);
    curr_pdg_actmaps = maps.PDG.cell_actmap_mat(:, :, :, currCell);
    box_size = size(curr_natmov_avgprojs, 1);

    %display all figures for current cell
    clf(fig)
    subplot('Position',[0.017 0.31 0.35 0.635])
    imagesc(curr_full_avgproj); colormap gray; axis square;
    set(gca,'XTick',[], 'YTick', [])
    circle(centers(currCell, 1), centers(currCell, 2), box_size/1.5);
    text(centers(currCell, 1)+21, centers(currCell, 2)-30, num2str(currCell), 'Color', 'w');
    temp_actmap = zeros(box_size, box_size, 3);
    for ii = 1:num_sessions
        subplot('Position', [0.3+(ii*0.085) 0.77 0.07 0.15])
        imagesc(curr_natmov_avgprojs(:, :, ii)); colormap gray; axis square;
        if ii == 1
            annotation('textbox', [0.37 0.94 0.07 0.01], 'String', 'NatMov average projections', 'FitBoxToText', 'on');
        end
        set(gca,'XTick',[], 'YTick', [])
    end
    for ii = 1:num_sessions
        subplot('Position', [0.3+(ii*0.085) 0.58 0.07 0.15])
        temp_actmap(:, :, 1) = curr_natmov_actmaps(:, :, ii);
        imagesc(temp_actmap); colormap gray; axis square;
        if ii == 1
            annotation('textbox', [0.37 0.75 0.07 0.01], 'String', 'NatMov activity maps', 'FitBoxToText', 'on');
        end
        set(gca,'XTick',[], 'YTick', [])
    end
    for ii = 1:num_sessions
        subplot('Position', [0.3+(ii*0.085) 0.39 0.07 0.15])
        imagesc(curr_pdg_avgprojs(:, :, ii)); colormap gray; axis square;
        if ii == 1
            annotation('textbox', [0.37 0.56 0.07 0.01], 'String', 'PDG average projections', 'FitBoxToText', 'on');
        end
        set(gca,'XTick',[], 'YTick', [])
    end
    for ii = 1:num_sessions
        subplot('Position', [0.3+(ii*0.085) 0.20 0.07 0.15])
        temp_actmap(:, :, 1) = curr_pdg_actmaps(:, :, ii);
        imagesc(temp_actmap); colormap gray; axis square;
        if ii == 1
            annotation('textbox', [0.37 0.37 0.07 0.01], 'String', 'PDG activity maps', 'FitBoxToText', 'on');
        end
        set(gca,'XTick',[], 'YTick', [])
        xlabel(sprintf('D%d', (ii-1)*7));
    end

    %slider
    UIelements.cellSlider = uicontrol(fig,'Style','slider','Units','normal',...
        'Position',[0.13 0.27 0.1 0.025],...
        'Min',1,'Max',numCells,...
        'Value',currCell,'SliderStep',[1/(numCells-1) 1/(numCells-1)],...
        'Parent',fig,...
        'Callback',@sliderCallback);
    
    %cell quality text indicators
    UIelements.qualityLabel = uicontrol(fig,'Style','text',...
        'String','Overall cell quality','Units','characters','Position',[57 17 20 1]);

    UIelements.qualitySelection = uicontrol(fig,'Style','text',...
        'String','0', 'Units','characters','Position',[67 8 20 3]);

    %Cell presence text indicators
    for yy = 1:num_sessions
        currfield = sprintf('presentS%d', yy);
        UIelements.(currfield) = uicontrol(fig, 'Style', 'text',...
        'String', '1', 'Units', 'characters', 'Position', [148+(yy-1)*31 10 10 1]);
    end

    %Cell quality pushbuttons
    UIelements.qualityPB5 = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normal', ...
        'Position', [0.13 0.20 0.05 0.025], 'Parent', fig,'String', '5',...
        'Callback', @QPB5Callback);

    UIelements.qualityPB4 = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normal', ...
        'Position', [0.13 0.16 0.05 0.025], 'Parent', fig,'String', '4',...
        'Callback', @QPB4Callback);

    UIelements.qualityPB3 = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normal', ...
        'Position', [0.13 0.12 0.05 0.025], 'Parent', fig,'String', '3',...
        'Callback', @QPB3Callback);

   UIelements. qualityPB2 = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normal', ...
        'Position', [0.13 0.08 0.05 0.025], 'Parent', fig,'String', '2',...
        'Callback', @QPB2Callback);

    UIelements.qualityPB1 = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normal', ...
        'Position', [0.13 0.04 0.05 0.025], 'Parent', fig,'String', '1',...
        'Callback', @QPB1Callback);

    %save pushbutton
    UIelements.savePB = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normal',...
        'Position', [0.27 0.16 0.06 0.05], 'Parent', fig, 'String', 'Save progress',...
        'Callback', @savePBCallback);
    
    %finished pushbutton
    UIelements.finishPB = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normal', ...
        'Position', [0.27 0.08 0.06 0.05], 'Parent', fig, 'String', 'Finished',...
        'Callback', @finishPBCallback);

    %cell precense pushbuttons
    for uu = 1:num_sessions
        currfield = sprintf('presentPB%d', uu);
        UIelements.(currfield) = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normal',...
        'Position', [0.395+(uu-1)*0.0855 0.08 0.05 0.04], 'Parent', fig, 'String', 'Present',...
        'Callback', {@presentPBCallback, uu});
    end

    function circle(x,y,r)                         
        hold on
        th = 0:pi/50:2*pi;
        xunit = r*cos(th) + x;
        yunit = r*sin(th) + y;
        plot(xunit, yunit, 'LineWidth', 2);
        hold off
    end

%Callbacks
    function sliderCallback(hObject,eventdata,handles)
        % change cell and rerun
        nextCell = uint16(get(hObject,'Value'));
        cellJudging_function(maps,nextCell,cellInfo,fig,stopFlag);
    end 
    %% Quality pushbutton callbacks
    function QPB5Callback(hObject, eventdata, handles)
        set(UIelements.qualitySelection, 'String', '5');
        cellInfo.quality(currCell) = 5;
    end

    function QPB4Callback(src, evt)
        set(UIelements.qualitySelection, 'String', '4');
        cellInfo.quality(currCell) = 4;
    end

    function QPB3Callback(src, evt)
        set(UIelements.qualitySelection, 'String', '3');
        cellInfo.quality(currCell) = 3;
    end

    function QPB2Callback(src, evt)
        set(UIelements.qualitySelection, 'String', '2');
        cellInfo.quality(currCell) = 2;
    end

    function QPB1Callback(src, evt)
        set(UIelements.qualitySelection, 'String', '1');
        cellInfo.quality(currCell) = 1;
    end
    %% presence pushbutton callback
    function presentPBCallback(src, evt, session)
        currSession = sprintf('presentS%d', session);
        currButton = UIelements.(currSession);
        currVal = get(currButton, 'String');
        switch currVal
            case '1'
                set(UIelements.(currSession), 'String', '0');
                cellInfo.presence(currCell, session) = 0;
            case '0'
                set(UIelements.(currSession), 'String', '1');
                cellInfo.presence(currCell, session) = 1;
        end
    end
    
%%
    function savePBCallback(src, evt)
        fprintf('Saving data...\n');
        save cellInfo cellInfo
    end

    function finishPBCallback(hObject, eventdata, handles)
        fprintf('Finished. Saving data...\n');
        save cellInfo cellInfo
        nextCell = uint16(get(hObject,'Value'))-1;
        stopFlag = 1;
        cellJudging_function(maps,nextCell,cellInfo,fig,stopFlag);
    end
end

