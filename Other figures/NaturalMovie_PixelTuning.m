function colorimg = NaturalMovie_PixelTuning(method, session, compress_flag)
%%% Written TM 201011, basd on Orientation pixel tuning code

if nargin == 0
    method = 'Full spectrum';         % 'Full spectrum'
    session = 1;
    compress_flag = 0;
end
% load file

disp('Select mat file: ')
[filename,~] = uigetfile('.mat');
load(filename)
tif_name = data.filename;
frameRate = data.frameRate;

% params
offTime = 5;                                         % time before stimulus (to analyze)
onTime = 30;                                          % grating time
repeats = 30;                                        % number of repeats

% calcualte additional params
repTime = offTime+onTime;
offFrames = round(offTime*frameRate);
onFrames = round(onTime*frameRate);
frameIdx = 1:onFrames; 

% average frames
fullMatrix = zeros(data.yPixels,data.xPixels,length(frameIdx),repeats);
for rep = 1:repeats
    progressbar(rep/repeats);
    stimFrame = ceil(((rep-1)*repTime+offTime)*data.frameRate);
    offIndex = stimFrame-offFrames+1:stimFrame;
    offResp = zeros(data.yPixels,data.xPixels,offFrames);
    for f = 1:offFrames
        offResp(:,:,f) = double(imread(tif_name,offIndex(f)));
    end
    onIndex = stimFrame+1:stimFrame+onFrames;
    onResp = zeros(data.yPixels,data.xPixels,onFrames);
    for f = 1:onFrames
        onResp(:,:,f) = double(imread(tif_name,onIndex(f)));
    end
    Resp = onResp-mean(offResp,3);
    Resp(Resp<0) = 0;

    fullMatrix(:,:,:,rep) = Resp;
end

meanResp = mean(fullMatrix, 4);
progressbar(1);
close all

fprintf('Calculating reliability...\n');
reliability = zeros(data.yPixels, data.xPixels);
for ii = 1:data.yPixels
    for jj = 1:data.xPixels
        jj
        curr_pixel = squeeze(fullMatrix(ii, jj, :, :));
        iterations = 10;
        CC_data_iterated = zeros(1, iterations);
        for qq = 1:iterations
            trialshuffle = randperm(repeats);
            trialselect1 = trialshuffle(1:floor(repeats/2));
            trialselect2 = trialshuffle(floor(repeats/2)+1:end);
            randavg1 = mean(curr_pixel(:, trialselect1), 2);
            randavg2 = mean(curr_pixel(:, trialselect2), 2);
            CC_data_iterated(qq) = corr(randavg1, randavg2);
        end
        reliability(ii, jj) = nanmean(CC_data_iterated);
    end
end


switch method
    case 'Full spectrum'
        if compress_flag == 1
            red_index = [251:300 1:112];
            red_gain = [0.5:(0.5/74):1 1:-(0.5/86):0.5];
            
            green_index = 76:175;
            green_gain = [0.5:(0.5/49):1 1:-(0.5/49):0.5];
            
            blue_index = 139:300;
            blue_gain = [0.5:(0.5/86):1 1:-(0.5/74):0.5];
        else
            
            red_index = [251:300 1:150];
            red_gain = [0:(1/99):1 1:-(1/99):0];
            
            green_index = 51:250;
            green_gain = [0:(1/99):1 1:-(1/99):0];
            
            blue_index = [151:300 1:50];
            blue_gain = [0:(1/99):1 1:-(1/99):0];
            
        end
        
        red_matrix = zeros(data.yPixels,data.xPixels);
        for i = 1:length(red_index)
            red_matrix = red_matrix+meanResp(:,:,red_index(i))*red_gain(i);
        end
        green_matrix = zeros(data.yPixels,data.xPixels);
        for i = 1:length(green_index)
            green_matrix = green_matrix+fullMatrix(:,:,green_index(i))*green_gain(i);
        end
        blue_matrix = zeros(data.yPixels,data.xPixels);
        for i = 1:length(blue_index)
            blue_matrix = blue_matrix+fullMatrix(:,:,blue_index(i))*blue_gain(i);
        end
        
        % normalize
        min_matrix = min([min(min(red_matrix)) min(min(green_matrix)) min(min(blue_matrix))]);
        red_matrix = red_matrix-min_matrix;
        green_matrix = green_matrix-min_matrix;
        blue_matrix = blue_matrix-min_matrix;
        max_matrix = max([max(max(red_matrix)) max(max(green_matrix)) max(max(blue_matrix))]);
        red_matrix = red_matrix/max_matrix;
        green_matrix = green_matrix/max_matrix;
        blue_matrix = blue_matrix/max_matrix;
        
        % concatenate
        colorimg = cat(3,red_matrix,green_matrix,blue_matrix);
    case 'Max'                          % NOT RECOMMENDED, PRODUCES COLORS THAT ARENT IN JET COLORMAP (DUE TO TAKING THE MEDIAN, IN ROIS WITH NOISY PIXEL COLORS)
        JCS003 = ResponseVisualizer;
        JCS003.importData;
        JCS003.setQualityThreshold(3);
        JCS003.cellSelection;
        qual = JCS003.getUse_cells;
        sessions = JCS003.getUse_sessions;
        good_cells = qual & sessions';
        

        [~, maxResp] = max(meanResp, [], 3);
        colorimg = ind2rgb(maxResp, jet(256));
        
        cellMasks = data.cellMasks;
        responsive = JCS003.RoiINFO.NatMov_Responsive_shuffle & JCS003.RoiINFO.NatMov_Responsive_shuffle;

        fullmask = zeros(data.yPixels, data.xPixels);
        for ii = 1:length(cellMasks)
            if good_cells(session, ii) && responsive(ii)
                fullmask = fullmask + poly2mask(cellMasks{ii}(:, 1),cellMasks{ii}(:, 2), data.yPixels, data.xPixels);
            end
        end

        for yy = 1:size(colorimg, 1)
            for xx = 1:size(colorimg, 2)
                if fullmask(yy, xx) ~= 1
                    colorimg(yy, xx, :) = 0; 
                end
            end
        end
        red = colorimg(:, :, 1);
        green = colorimg(:, :, 2);
        blue = colorimg(:, :, 3);
        for ii = 1:length(cellMasks)
            currmask = poly2mask(cellMasks{ii}(:, 1),cellMasks{ii}(:, 2), data.yPixels, data.xPixels);
            redvalues = red(currmask);
            greenvalues = green(currmask);
            bluevalues = blue(currmask);
            hexvalues = cell(1, length(redvalues));
            for pp = 1:length(redvalues)
                curr_RGB = [redvalues(pp) greenvalues(pp) bluevalues(pp)];
                hexvalues{pp} = rgb2hex(curr_RGB);
            end
            uniquevalues = unique(hexvalues);
            num_pixels = zeros(1, length(uniquevalues));
            for uu = 1:length(uniquevalues)
                num_pixels(uu) = sum(strcmp(hexvalues, uniquevalues(uu)));
            end
            
            [~, mostcommon_idx] = max(num_pixels);
            mostcommon_hex = uniquevalues(mostcommon_idx);           
            mostcommon_RGB = hex2rgb(mostcommon_hex{1});
            
            red(currmask) = mostcommon_RGB(1);
            green(currmask) = mostcommon_RGB(2);
            blue(currmask) = mostcommon_RGB(3);       
        end
        colorimg(:, :, 1) = red;
        colorimg(:, :, 2) = green;
        colorimg(:, :, 3) = blue;
end




% addmasks = 1;
% colorimg = data.colorimg;

% display and save figure
if strcmp(method, 'Full spectrum')
    gain = 3;
    gain_mat = reliability*gain;
    for i = 1:3
        colorimg(:, :, i) = colorimg(:, :, i).*gain_mat;
    end

%     JCS003 = ResponseVisualizer;
%     JCS003.importData;
%     JCS003.setQualityThreshold(3);
%     JCS003.cellSelection;
%     qual = JCS003.getUse_cells;
%     sessions = JCS003.getUse_sessions;
%     good_cells = qual & sessions';
%     
%     
%     cellMasks = data.cellMasks;
%     responsive = JCS003.RoiINFO.NatMov_Responsive_shuffle & JCS003.RoiINFO.NatMov_Responsive_thresh;
%     fullmask = zeros(data.yPixels, data.xPixels);
%     for ii = 1:length(cellMasks)
%         if good_cells(session, ii) && responsive(ii)
%             fullmask = fullmask + poly2mask(cellMasks{ii}(:, 1),cellMasks{ii}(:, 2), data.yPixels, data.xPixels);
%         end
%     end 
%     
%     if addmasks == 1
%         new_Masks = importdata('new_cellMask.mat');
%         for ll = 1:length(new_Masks)
%             fullmask = fullmask + poly2mask(new_Masks{ll}(:, 1), new_Masks{ll}(:, 2), data.yPixels, data.xPixels);
%         end
%     end
    
%     roi_display_gain = [3 3 3];
%     nonroi_display_gain = [1 1 1];
%     fullmask = logical(fullmask);
%     for i = 1:3
%         currplane = colorimg(:, :, i);
%         currplane(fullmask) = currplane(fullmask)*roi_display_gain(i);
%         currplane(~fullmask) = currplane(~fullmask)*nonroi_display_gain(i);
%         colorimg(:,:,i) = currplane;
%     end
%     colorimg(colorimg>1) = 1;
end



% save matrix
data.colorimg = colorimg;
save(filename,'data');


figure
image(colorimg)
axis square
set(gcf,'color',[1 1 1])
set(gcf,'Position',[300 100 1000 1000])
box off
set(gca, 'xtick', []);
set(gca, 'ytick', []);
% saveas(gcf,'mergeOri')
% shg
% pause
% close
end

