function [layerVec, boundary234_Pos, boundary45_Pos, boundary56_Pos, window_midpoints_out, ROIdensity_out] = drawLayerBoundaries_V2(um_per_pixel, method)
%%% Allow user to draw layer boundaries using activity map
%%% Assign layer to each ROI based on centroid
%%% Written TM MG last updated 190422

if nargin < 2
    method = 'auto';
end

% load file
fprintf('Load dummy data file (exact copy of real data file, NAMED DIFFERENTLY) for rotation and boundary drawing... \n...New ROIs should be drawn purely using automatic detection on this dummy datafile to draw layer boundaries in unbiased fashion\n');
[filename,pathname] = uigetfile('*.mat','Load data');
cd(pathname)
dummy_data = importdata(filename);

%Load in datafile with real ROIs for the rest
fprintf('Load real data file containing original cell masks\n')
[filename_real,pathname_real] = uigetfile('*.mat','Load data');
cd(pathname)
real_data = importdata(filename_real);



Flag234=0;
Flag45=0;
Flag56=0;

if strcmp(method, 'auto')        %attempt to autogenerate layer boundaries

    %find vertical axis and calculate angle for rotation
    figure
    imagesc(dummy_data.activity_map)
    axis square
    colormap gray
    set(gcf,'position',[400 150 1000 850])
    fprintf('Draw vertical axis of window (cortical column)\n');
    title('Draw vertical axis (cortical column)')
    V = drawline;
    pause
    vertical = V.Position;
    close
    if vertical(1,2) > vertical(2,2)            %sort by ascending y coord
        vertical = [vertical(2,:); vertical(1,:)];
    end
    y1 = vertical(1,2);
    y2 = vertical(2,2);
    x1 = vertical(1,1);
    x2 = vertical(2,1);
    theta = atan((y1 - y2)/(x1 - x2))*180/pi;
    
    %redefine ROIs
    B_DefineROI();
    
    dummy_data = importdata(filename);
    
    maskImg = false(dummy_data.yPixels,dummy_data.xPixels);

    for ii = 1:length(dummy_data.cellMasks)
        position = dummy_data.cellMasks{ii};
        mask = poly2mask(position(:,1),position(:,2),dummy_data.yPixels,dummy_data.xPixels);

        maskImg(mask) = true;

    end

    %rotate avg proj and activity map, save to new data file and define new ROIs on this file to be used for ROI analysis 
    rotated_avg = imrotate(dummy_data.avg_projection, theta, 'bilinear', 'loose');      %edit these to rotate by the apropriate angle based on the drawn vertical
    rotated_act = imrotate(dummy_data.activity_map, theta, 'bilinear', 'loose');
    rotated_cellMasks = imrotate(maskImg, theta, 'bilinear', 'loose');
    dummy_data.cellMasks = bwboundaries(transpose(rotated_cellMasks), 4)';
    dummy_data.avg_projection = rotated_avg;
    dummy_data.activity_map = rotated_act;
    dummy_data.xPixels = size(rotated_avg, 1);
    dummy_data.yPixels = size(rotated_avg, 2);
    
    %Calculate centroids of all ROIs just created
    centroids = zeros(length(dummy_data.cellMasks), 2);
    for cell = 1:length(dummy_data.cellMasks)
        centroids(cell, :) = mean(dummy_data.cellMasks{cell});
    end
    
    %iterate through the cortex in certtain size windows and find average ROI size for all centroids that fall within each window
    window_size =  65;       %number of pixels per window
    num_blocks = ceil(dummy_data.xPixels/window_size);
    avg_ROIsize = zeros(1, num_blocks);
    
    for n = 1:num_blocks
        bound1 = (n-1)*window_size+1;
        bound2 = n*window_size;
        
        cells_in_bounds = (centroids(:, 1) > bound1) & (centroids(:, 1) < bound2);
        cells_to_avg = dummy_data.cellMasks(cells_in_bounds);
        
        roiSizes = zeros(1, length(cells_to_avg));
        
        for cell = 1:length(cells_to_avg)
            roiSizes(cell) = length(cells_to_avg{cell});
        end
        
        avg_ROIsize(n) = mean(roiSizes);
    end
    
    
    %iterate through again in (different?) size windows to determine density of cells in each window
    ROIdensity = zeros(1, num_blocks);
    
    for n = 1:num_blocks
        bound1 = (n-1)*window_size+1;
        bound2 = n*window_size;
        
        cells_in_bounds = (centroids(:, 1) > bound1) & (centroids(: ,1) < bound2);
        ROIdensity(n) = sum(cells_in_bounds)/(window_size*dummy_data.yPixels);
    end
    
    figure; 
    yyaxis left
    plot((1:num_blocks)*window_size - floor(window_size/2), avg_ROIsize, 'LineWidth', 1.5);
    ylabel('Average ROI size (pixels)')
    ylim([50 80])
    yyaxis right
    plot((1:num_blocks)*window_size - floor(window_size/2), ROIdensity, 'LineWidth', 1.5); 
    ylabel('Average ROI density (1/pixel^2)')
    ylim([0 0.0025])
    xlabel('Pixels across vertical axis (approximate)')
    xlim([0 850])
    axis square
    
    window_midpoints_out = (1:num_blocks)*window_size - floor(window_size/2);
    ROIdensity_out = ROIdensity;
    
dummy_w = dummy_data.xPixels;
real_w = real_data.xPixels;
difference_w = (dummy_w-real_w)/2;

dummy_h = dummy_data.yPixels;
real_h = real_data.yPixels;
difference_h = (dummy_h-real_h)/2;
    
%Determine layer2/3-4 boundary based on above characteristics

cf = 0;

L4um = 140;     % layer thickness in microns
L5um = 150;

L4thickness = (L4um/um_per_pixel)*(dummy_data.xPixels/real_data.xPixels);            % pixel thickness of layer 4, assume layer 4 is about 140um thick, scale by increase in pixel size due to rotation
L5thickness = (L5um/um_per_pixel)*(dummy_data.xPixels/real_data.xPixels);


[~, L4window] = max(ROIdensity);                    %using only ROI density curve because it seems to be most reliable and sufficient
L4bound = round((L4window+cf)*window_size - L4thickness/2);
    
    %create new boundary line and rotate it to original map
    Lv = zeros(2);
    Lv(1, 1) = L4bound;                             
    Lv(1, 2) = L4bound;
    Lv(2, 1) = find(rotated_act(:, L4bound), 1)+1;
    Lv(2, 2) = find(rotated_act(:, L4bound), 1, 'last')-1;
    
    Lc = zeros(2,1);
    Lc(1, 1) = L4bound;
    Lc(2, 1) = mean([Lv(2,1) Lv(2,2)]);
    Lc = repmat(Lc, 1, 2);
    R = [cos(theta*pi/180) -sin(theta*pi/180); sin(theta*pi/180) cos(theta*pi/180)];
    vo = R*(Lv - Lc) + Lc;

    boundary234_Pos(:,1) = vo(1, :);
    boundary234_Pos(:,2) = vo(2, :);
    
    % correct for width and height changes due to rotation
    boundary234_Pos(:, 1) = boundary234_Pos(:, 1) - difference_w;
    boundary234_Pos(:, 2) = boundary234_Pos(:, 2) - difference_h;
    
    Flag234 = 1;
    
    
    %Determine layer4-5 boundary based on above characteristics
L5bound = round((L4window+cf)*window_size + L4thickness/2);

    Lv = zeros(2);
    Lv(1, 1) = L5bound;                             
    Lv(1, 2) = L5bound;
    Lv(2, 1) = find(rotated_act(:, L5bound), 1)+1;
    Lv(2, 2) = find(rotated_act(:, L5bound), 1, 'last')-1;
    
    Lc = zeros(2,1);
    Lc(1, 1) = L5bound;
    Lc(2, 1) = mean([Lv(2,1) Lv(2,2)]);
    Lc = repmat(Lc, 1, 2);
    R = [cos(theta*pi/180) -sin(theta*pi/180); sin(theta*pi/180) cos(theta*pi/180)];
    vo = R*(Lv - Lc) + Lc;

    boundary45_Pos(:,1) = vo(1, :);
    boundary45_Pos(:,2) = vo(2, :);
    
    boundary45_Pos(:, 1) = boundary45_Pos(:, 1) - difference_w;
    boundary45_Pos(:, 2) = boundary45_Pos(:, 2) - difference_h;
%     
    Flag45 = 1;
    
L6bound = round(L5bound + L5thickness);

    Lv = zeros(2);
    Lv(1, 1) = L6bound;                             
    Lv(1, 2) = L6bound;
    Lv(2, 1) = find(rotated_act(:, L6bound), 1)+1;
    Lv(2, 2) = find(rotated_act(:, L6bound), 1, 'last')-1;
    
    Lc = zeros(2,1);
    Lc(1, 1) = L6bound;
    Lc(2, 1) = mean([Lv(2,1) Lv(2,2)]);
    Lc = repmat(Lc, 1, 2);
    R = [cos(theta*pi/180) -sin(theta*pi/180); sin(theta*pi/180) cos(theta*pi/180)];
    vo = R*(Lv - Lc) + Lc;

    boundary56_Pos(:,1) = vo(1, :);
    boundary56_Pos(:,2) = vo(2, :);
    
    boundary56_Pos(:, 1) = boundary56_Pos(:, 1) - difference_w;
    boundary56_Pos(:, 2) = boundary56_Pos(:, 2) - difference_h;
%     
    Flag56 = 1;
end



LoopFlag = 0;

while LoopFlag == 0
    fprintf('Manual layer boundary adjustment (if necessary). Select boundary to redraw.\n');
    [BoundaryType, ~] = listdlg('ListString', {'L2/3-4', 'L4-5', 'L5-6', 'Done'});
    
    figure
    imagesc(real_data.activity_map)
    if Flag234 == 1
        hold on
        plot(boundary234_Pos(:,1), boundary234_Pos(:,2), 'r-');
    end
    if Flag45 == 1
        hold on
        plot(boundary45_Pos(:,1), boundary45_Pos(:,2), 'r-');
    end
     if Flag56 == 1
        hold on
        plot(boundary56_Pos(:,1), boundary56_Pos(:,2), 'r-');
    end
    axis square
    colormap('gray')
    set(gcf,'position',[400 150 1000 850])
    
    switch BoundaryType
        case 1
            fprintf('Draw boundary between L2/3 and L4, press any key when done\n')
            title('Draw boundary between L2/3 and L4, press any key when done')
            H = drawline;
            pause
            boundary234_Pos = H.Position;
            close
            Flag234 = 1;
        case 2
            fprintf('Draw boundary between L4 and L5, press any key when done\n')
            title('Draw boundary between L4 and L5, press any key when done')
            H = drawline;
            pause
            boundary45_Pos = H.Position;
            close
            Flag45 = 1;
        case 3
            fprintf('Draw boundary between L5 and L6, press any key when done\n')
            title('Draw boundary between L5 and L6, press any key when done')
            H = drawline;
            pause
            boundary56_Pos = H.Position;
            close
            Flag56 = 1;
        case 4
            LoopFlag = 1;
    end
end
            

% Arrange by ascending Y-coordinate
if Flag234 == 1 && (boundary234_Pos(1,2) > boundary234_Pos(2,2))
    boundary234_Pos = [boundary234_Pos(2,:); boundary234_Pos(1,:)];
end
if Flag45 == 1 && (boundary45_Pos(1,2) > boundary45_Pos(2,2))
    boundary45_Pos = [boundary45_Pos(2,:); boundary45_Pos(1,:)];
end
if Flag56 == 1 && (boundary56_Pos(1,2) > boundary56_Pos(2,2))
    boundary56_Pos = [boundary56_Pos(2,:); boundary56_Pos(1,:)];
end



if Flag234 == Flag45 == Flag56 == 1
% Assign layers
    for cell = 1:length(real_data.cellMasks)

        % Extract centroid
        currROI = real_data.cellMasks{cell};
        currCentroid = mean(currROI);
        xc = currCentroid(1);
        yc = currCentroid(2);

        % Determine if left of boundary 3
        x1 = boundary56_Pos(1,1);
        x2 = boundary56_Pos(2,1);
        y1 = boundary56_Pos(1,2);
        y2 = boundary56_Pos(2,2);
        distCentroid = (xc-x1)*(y2-y1)-(yc-y1)*(x2-x1);

        if distCentroid > 0
            layerVec{cell} = 'L6';
        else
            % Determine if left of boundary 2
            x1 = boundary45_Pos(1,1);
            x2 = boundary45_Pos(2,1);
            y1 = boundary45_Pos(1,2);
            y2 = boundary45_Pos(2,2);
            distCentroid = (xc-x1)*(y2-y1)-(yc-y1)*(x2-x1);
            if distCentroid > 0
                layerVec{cell} = 'L5';
            else
                x1 = boundary234_Pos(1,1);
                x2 = boundary234_Pos(2,1);
                y1 = boundary234_Pos(1,2);
                y2 = boundary234_Pos(2,2);
                distCentroid = (xc-x1)*(y2-y1)-(yc-y1)*(x2-x1);
                if distCentroid > 0
                    layerVec{cell} = 'L4';
                else layerVec{cell} = 'L2/3';
                end
            end
        end
    end
end


if Flag234 == Flag45 == 1 && Flag56 == 0 
% Assign layers
    for cell = 1:length(real_data.cellMasks)

        % Extract centroid
        currROI = real_data.cellMasks{cell};
        currCentroid = mean(currROI);
        xc = currCentroid(1);
        yc = currCentroid(2);

        % Determine if left of boundary 3
        x1 = boundary45_Pos(1,1);
        x2 = boundary45_Pos(2,1);
        y1 = boundary45_Pos(1,2);
        y2 = boundary45_Pos(2,2);
        distCentroid = (xc-x1)*(y2-y1)-(yc-y1)*(x2-x1);

        if distCentroid > 0
            layerVec{cell} = 'L5';
        else
            % Determine if left of boundary 2
            x1 = boundary234_Pos(1,1);
            x2 = boundary234_Pos(2,1);
            y1 = boundary234_Pos(1,2);
            y2 = boundary234_Pos(2,2);
            distCentroid = (xc-x1)*(y2-y1)-(yc-y1)*(x2-x1);
            if distCentroid > 0
                layerVec{cell} = 'L4';
            else layerVec{cell} = 'L2/3';
            end
        end
    end
end

if Flag45 == Flag56 == 1 && Flag234 == 0
    
% Assign layers
    for cell = 1:length(real_data.cellMasks)

        % Extract centroid
        currROI = real_data.cellMasks{cell};
        currCentroid = mean(currROI);
        xc = currCentroid(1);
        yc = currCentroid(2);

        % Determine if left of boundary 3
        x1 = boundary56_Pos(1,1);
        x2 = boundary56_Pos(2,1);
        y1 = boundary56_Pos(1,2);
        y2 = boundary56_Pos(2,2);
        distCentroid = (xc-x1)*(y2-y1)-(yc-y1)*(x2-x1);

        if distCentroid > 0
            layerVec{cell} = 'L6';
        else
            % Determine if left of boundary 2
            x1 = boundary45_Pos(1,1);
            x2 = boundary45_Pos(2,1);
            y1 = boundary45_Pos(1,2);
            y2 = boundary45_Pos(2,2);
            distCentroid = (xc-x1)*(y2-y1)-(yc-y1)*(x2-x1);
            if distCentroid > 0
                layerVec{cell} = 'L5';
            else layerVec{cell} = 'L4';
            end
        end
    end
end

figure
imagesc(real_data.avg_projection);
hold on
colormap gray; axis square;
for cell = 1:length(real_data.cellMasks)
    currROI = real_data.cellMasks{cell};
    if strcmp(layerVec{cell}, 'L2/3')
        color = [1 0 0];
    elseif strcmp(layerVec{cell}, 'L4')
        color = [0 1 0];
    elseif strcmp(layerVec{cell}, 'L5')
        color = [0 0 1];
    else
        color = [0 0 0];
    end
    pdata = patch(currROI(:,1),currROI(:,2), color);
    pdata.EdgeAlpha = 0;
    if ~strcmp(layerVec{cell}, 'L6')
        pdata.FaceAlpha = 0.5;
    else
        pdata.FaceAlpha = 0;
    end
end


end

        
            
    
    
    
