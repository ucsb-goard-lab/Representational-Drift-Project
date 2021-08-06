function ReceptiveField_preProcessing(fn_Stimdat, pn_Stimdat, fn_data, pn_data, overwrite_flag)
%Step 1/2: Process each direction using this first, then we'll combine it later. 
% Preprocess and sort receptive field mapping response data
%Written KS
%Updated TM 190718

% Get the information from the Stimdata
if nargin == 0
    [fn_Stimdat,pn_Stimdat] = uigetfile('.mat');
    [fn_data,pn_data] = uigetfile('.mat');
end
if nargin < 5
    overwrite_flag = 'Yes';
end

Stimdata = importdata(fn_Stimdat);
data     = importdata(fn_data);

if isfield(data, 'isVisuallyResponsive') && strcmp(overwrite_flag, 'No')
    disp('File has already been pre-processed. Skipping...');
    return
end
if ~isfield(data,'spikes')
    data = spikeInference(fn_data,'Yes');
end

% Get the data that you need

zThresh = 2;

on_time = Stimdata.on_time;
off_time = Stimdata.off_time;
repeats = Stimdata.repeats;
bar_centers = Stimdata.centers;
num_locations = size(bar_centers,2);
fs = 10;

% Convert to frames
on_frames = on_time*fs;
off_frames = off_time*fs;

% Calculate additional times
rep_frames = off_frames + (on_frames * num_locations);

% Sort the data appropriately (reps and shuffled) (cells, on_frames, locs,reps)
RespVec = zeros(size(data.DFF,1),on_frames,num_locations,repeats);

for rep = 1:repeats
    curr_frame = (rep-1)*rep_frames;
    off_resp(:,:,rep) = data.spikes(:,curr_frame+1 : curr_frame+off_frames);
    for loc = 1:num_locations
        curr_frame = (rep-1)*rep_frames + (loc-1)*on_frames + off_frames;
        on_resp(:,:,loc,rep) = data.spikes(:,curr_frame+1:curr_frame+on_frames);
        RespVec(:,:,loc,rep) = on_resp(:,:,loc,rep) - mean(off_resp(:,:,rep),2);
    end
end

% Unscramble the data using Stimdat

for rep = 1:repeats
    [~,sorting_vector] = sort(bar_centers(rep,:));
    data.RespVec(:,:,:,rep) = RespVec(:,:,sorting_vector,rep);
    on_resp(:,:,:,rep) = on_resp(:,:,sorting_vector,rep);
end

% Check for visual responsiveness
off_dev = std(mean(off_resp,3),[],2); % mean across frames and rep
mean_on_resp = squeeze(mean(mean(on_resp,4),2)); % mean across frames and rep
max_on_resp = max(mean_on_resp,[],2); % preferred location resp
zScoreVec = max_on_resp ./ off_dev;
data.isVisuallyResponsive = zScoreVec > zThresh;

disp(['Percent visually responsive: ' num2str(mean(data.isVisuallyResponsive)*100) '%'])
save(fn_data,'data')

end