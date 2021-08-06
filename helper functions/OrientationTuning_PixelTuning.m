%%% OrientationTuning_PixelTuning
%%% Written MG 160613

% load file
clear
close all
disp('Select mat file: ')
[filename,~] = uigetfile('.mat');
load(filename)
tif_name = data.filename;
frameRate = data.frameRate;

% Gratings params
preTime = 2;                                         % recovery time (to ignore)
offTime = 2;                                         % time before stimulus (to analyze)
onTime = 2;                                          % grating time
numOri = 12;                                         % number of gratings
repeats = 8;                                        % number of repeats

% calcualte additional params
oriTime = preTime+offTime+onTime;
repTime = oriTime*numOri;
offFrames = round(offTime*frameRate);
onFrames = round(onTime*frameRate);
oriIndex = rem([1:numOri]-1,numOri/2)+1;  % collapse across opposite oris

% average frames
oriMatrix = zeros(data.yPixels,data.xPixels,numOri/2);
for rep = 1:repeats
    progressbar(rep/repeats);
    for ori = 1:numOri
        stimFrame = ceil(((rep-1)*repTime+(ori-1)*oriTime+preTime+offTime)*data.frameRate);
        offIndex = [stimFrame-offFrames:stimFrame-1];
        offResp = zeros(data.yPixels,data.xPixels,offFrames);
        for f = 1:offFrames
            offResp(:,:,f) = double(imread(tif_name,offIndex(f)));
        end
        onIndex = [stimFrame:stimFrame+onFrames-1];
        onResp = zeros(data.yPixels,data.xPixels,onFrames);
        for f = 1:onFrames
            onResp(:,:,f) = double(imread(tif_name,onIndex(f)));
        end
        oriResp = mean(onResp,3)-mean(offResp,3);
        oriResp(oriResp<0) = 0;
        
        oriMatrix(:,:,oriIndex(ori)) = oriMatrix(:,:,oriIndex(ori))+oriResp;
    end
end
progressbar(1);
close all

% make color matrices
red_index = [6 1 2];
red_gain = [0.5 1 0.5];
red_matrix = zeros(data.yPixels,data.xPixels);
for i = 1:length(red_index)
    red_matrix = red_matrix+oriMatrix(:,:,red_index(i))*red_gain(i);
end

green_index = [2 3 4];
green_gain = [0.5 1 0.5];
green_matrix = zeros(data.yPixels,data.xPixels);
for i = 1:length(green_index)
    green_matrix = green_matrix+oriMatrix(:,:,green_index(i))*green_gain(i);
end

blue_index = [4 5 6];
blue_gain = [0.5 1 0.5];
blue_matrix = zeros(data.yPixels,data.xPixels);
for i = 1:length(blue_index)
    blue_matrix = blue_matrix+oriMatrix(:,:,blue_index(i))*blue_gain(i);
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
merge_matrix = cat(3,red_matrix,green_matrix,blue_matrix);

% save matrix
data.merge_matrix = merge_matrix;
save(filename,'data');

% display and save figure
display_gain = [2 2 2];
for i = 1:3
    merge_matrix(:,:,i) = merge_matrix(:,:,i)*display_gain(i);
end
merge_matrix(merge_matrix>1) = 1;
image(merge_matrix)
axis square
set(gcf,'color',[1 1 1])
set(gcf,'Position',[300 100 1100 1000])
saveas(gcf,'mergeOri')
shg
pause
close

