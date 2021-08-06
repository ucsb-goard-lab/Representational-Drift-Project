
%% example cell maps
all = ChronicImaging_Mouse;
all.importData;                     % import datafiles for desired field ('RespData' & 'RoiINFO')
all.evaluateROIs;                   % run visual inspection preprocessor ('Visual_Inspection_A') and make cell maps available
num_sessions = all.num_sessions;
num_cells = all.num_cells;
neuron = 1;   


% run this section to flip the avg proj images
% fixing image artifacts by blowing up and rotating
box_size = size(all.Vis_cellInspection_maps.NatMov.cell_avgproj_mat, 1);
NatMov_maps = zeros(box_size, box_size, num_sessions, num_cells);                   % avg projections for natmov and pdg (act maps remain unaltered
PDG_maps = zeros(box_size, box_size, num_sessions, num_cells);
for ii = 1:num_cells
    ii
    for kk = 1:num_sessions
        currmap = all.Vis_cellInspection_maps.NatMov.cell_avgproj_mat(:, :, kk, ii);
        currmap = imresize(currmap, [1000 1000]);
        currmap = imrotate(currmap, -90);
        currmap = imresize(currmap, [36 36]);
        NatMov_maps(:, :, kk, ii) = currmap;
        
        currmap = all.Vis_cellInspection_maps.PDG.cell_avgproj_mat(:, :, kk, ii);
        currmap = imresize(currmap, [1000 1000]);
        currmap = imrotate(currmap, -90);
        currmap = imresize(currmap, [36 36]);
        PDG_maps(:, :, kk, ii) = currmap;
    end
end      

% run this section to keep them as they are straight out of the preprocessor
NatMov_maps = all.Vis_cellInspection_maps.NatMov.cell_avgproj_mat;
PDG_maps = all.Vis_cellInspection_maps.PDG.cell_avgproj_mat;

% making the figures
q = 5;
neuron = 295;
figure
for ii = neuron:all.num_cells % ii = cell;
    if all.RoiINFO.quality(ii) == q
        for kk = 1:num_sessions
            subplot(4, 7, kk)
            imagesc(NatMov_maps(:, :, kk, ii));
            colormap gray
            axis square
        end
        for kk = 1:num_sessions
            subplot(4, 7, kk+7)
            imagesc(imgaussfilt(all.Vis_cellInspection_maps.NatMov.cell_actmap_mat(:, :, kk, ii)));
            colormap gray
            axis square
        end
        for kk = 1:num_sessions
            subplot(4, 7, kk+14)
            imagesc(PDG_maps(:, :, kk, ii));
            colormap gray
            axis square
        end
        for kk = 1:num_sessions
            subplot(4, 7, kk+21)
            imagesc(imgaussfilt(all.Vis_cellInspection_maps.PDG.cell_actmap_mat(:, :, kk, ii)));
            colormap gray
            axis square
        end
        title(sprintf('cell %d, quality %d', ii, all.RoiINFO.quality(ii)));
        pause
    end
end


%% reliability on D0 for example field
JCS003 = StabilityAnalyzer;             % initialize Analyzer object
JCS003.importData;                      % import datafiles for example field ('RespData', 'RoiINFO', 'StabilityData')
JCS003.setQualityThreshold(3);          % set ROI quality threshold
JCS003.cellSelection;                   % set inclusion criteria (e.g. select 'quality' & 'sessions' for all neurons regardless of responsiveness)
qual = JCS003.getUse_cells;
sessions = JCS003.getUse_sessions;

% all quality neurons present on first session PDG vs MOV
figure
hold on
histogram(JCS003.StabilityData.PDG.CCs(1, qual' & sessions(:, 1)), 0:0.1:1, 'FaceColor', 'b');
histogram(JCS003.StabilityData.NatMov.CCs(1, qual' & sessions(:, 1)), 0:0.1:1, 'FaceColor', 'r');
ylim([0 120])

% PDG: all quality neurons present on first session & all quality, visually responsive neurons present on first session 
figure
hold on
histogram(JCS003.StabilityData.PDG.CCs(1, qual' & sessions(:, 1)), 0:0.1:1, 'FaceColor', 'b', 'FaceAlpha', 1);
histogram(JCS003.StabilityData.PDG.CCs(1, qual' & sessions(:, 1) & JCS003.RoiINFO.PDG_Responsive_thresh'), 0:0.1:1, 'FaceColor', 'c', 'FaceAlpha', 1);

% MOV: all quality neurons present on first session & all quality, visually responsive neurons present on first session 
figure
hold on
histogram(JCS003.StabilityData.NatMov.CCs(1, qual' & sessions(:, 1)), 0:0.1:1, 'FaceColor', 'r', 'FaceAlpha', 1);
histogram(JCS003.StabilityData.NatMov.CCs(1, qual' & sessions(:, 1) & JCS003.RoiINFO.NatMov_Responsive_thresh'), 0:0.1:1, 'FaceColor', 'm', 'FaceAlpha', 1);


%% example cell RDI curve
% utilizes Analyzer object created above
neuron = 86;         
NatMovcurve = JCS003.StabilityData.NatMov.RDI(:, neuron);
PDGcurve = JCS003.StabilityData.PDG.RDI(:, neuron);
figure
hold on
plot(PDGcurve, 'color', [0 0.7 0.7])
plot(NatMovcurve, 'color', [0.8 0.5 0.5])
refline(0, JCS003.StabilityData.PDG.RDI_control(neuron));
refline(0, JCS003.StabilityData.NatMov.RDI_control(neuron));
ylim([-0.2 1])
xlim([0 JCS003.num_sessions])
axis square; box off


%% Reliability and responsivity
all = StabilityAnalyzer;
all.importData;         % import datafiles for first mouse ('RespData', 'RoiINFO', 'StabilityData')
all.addData();          % add datafiles for subsequent mice (all consecutive sessions)
all.addData([1 1 0 1 1 1]);         % if necessary, use input arg to indicate weeks for which data is present, only necessary for mice with non-consecutive sessions. 
                                    % e.g. [1 1 1 0 1 1 1] for data corresponding to sessions 1 2 3 5 6 7, [1 1 0 1 1 1] for sessions 1 2 4 5 6, etc.
all.setQualityThreshold(3);
all.cellSelection;          % select 'quality' & 'presence' here
qual = all.getUse_cells;
sessions = all.getUse_sessions;

sessions(isnan(sessions)) = 0;            % nan values converted to zeros for logical vector purposes

good_cells = qual' & sessions;                % cells of sufficient quality & presence on any given session
D0_good_cells = good_cells(:, 1);                             % quality cells present on D0

all.cellSelection;               % select quality + all responsiveness criteria (i.e. PDG_Responsive_thresh & NatMov_Responsive_thresh)
dual_responsive = all.getUse_cells';       %cells responsive to both stimuli

% generating vector of cells per mouse
num_session_allmice = [7 7 5 6 6 5 7 7 7 7 5 7 7];                %number of sessions for each mouse, in order in which they were added to the Analyzer
num_cells_allmice = all.getCellList;

ct = 1;
for ii = 1:length(num_session_allmice)
    cellstart_idx(ii) = ct;                                    % starting index for the first cell in every mouse
    ct = ct + num_cells_allmice(ii);
end

% %%%PDG vs NatMov reliability on D0, using good cells present on D0 and responsive to both stimuli
% figure
% scatter(JCS003.StabilityData.PDG.CCs(1, D0_good_cells & dual_responsive), JCS003.StabilityData.NatMov.CCs(1, D0_good_cells & dual_responsive), 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4);
% [r, p] = corr(JCS003.StabilityData.PDG.CCs(1, D0_good_cells & dual_responsive)', JCS003.StabilityData.NatMov.CCs(1, D0_good_cells & dual_responsive)');
% title(sprintf('Reliability (dual responsive), r = %.2f, p = %.4f', r, p));
% xlabel('PDG')
% ylabel('NatMov')
% box off
% axis square
% refline(1, 0)

%%%PDG vs NatMov reliability on D0, using good cells present on D0 (all neurons regardless of responsiveness)
PDG_Responsive = all.RoiINFO.PDG_Responsive_thresh;  
NatMov_Responsive = all.RoiINFO.NatMov_Responsive_thresh;    
PDG_only = qual & PDG_Responsive & ~NatMov_Responsive;
NatMov_only = qual & ~PDG_Responsive & NatMov_Responsive;
both = qual & PDG_Responsive & NatMov_Responsive;
none = qual & ~PDG_Responsive & ~NatMov_Responsive;

% colored by responsiveness to stimuli
figure
hold on
scatter(all.StabilityData.PDG.CCs(1, qual & none), all.StabilityData.NatMov.CCs(1, qual & none), 30, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerFaceColor', [0.7 0.7 0.7]);
scatter(all.StabilityData.PDG.CCs(1, qual & PDG_only), all.StabilityData.NatMov.CCs(1, qual & PDG_only), 30, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerFaceColor', [0 0.7 0.7]);
scatter(all.StabilityData.PDG.CCs(1, qual & NatMov_only), all.StabilityData.NatMov.CCs(1, qual & NatMov_only), 30, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerFaceColor', [0.8 0.5 0.5]);
scatter(all.StabilityData.PDG.CCs(1, qual & both), all.StabilityData.NatMov.CCs(1, qual & both), 30, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [1.0 0.9 0]);

[r, p] = corr(all.StabilityData.PDG.CCs(1, qual)', all.StabilityData.NatMov.CCs(1, qual)');
title(sprintf('Reliability (all), r = %.2f, p = %.4f', r, p));
xlabel('PDG')
ylabel('NatMov')
box off
axis square
refline(1, 0)
xlim([0 1])
ylim([0 1])

%% RDI curves, all neurons
%uses StabilityAnalyzer object generated above using all imaging fields
all = StabilityAnalyzer;
all.importData;         % import datafiles for first mouse ('RespData', 'RoiINFO', 'StabilityData')
all.addData();          % add datafiles for subsequent mice (all consecutive sessions)
all.addData([1 1 0 1 1 1]);         % if necessary, use input arg to indicate weeks for which data is present, only necessary for mice with non-consecutive sessions. 
                                    % e.g. [1 1 1 0 1 1 1] for data corresponding to sessions 1 2 3 5 6 7, [1 1 0 1 1 1] for sessions 1 2 4 5 6, etc.
all.setQualityThreshold(3);
all.cellSelection;
qual = all.getUse_cells;
sessions = all.getUse_sessions;

sessions(isnan(sessions)) = 1;            % nan values converted to ones for logical vector purposes
present = sum(sessions, 2) == all.num_sessions;
good_cells = qual' & present;

figure
all.RDIplotter('PDG', 'Average');
hold on
all.RDIplotter('NatMov', 'Average');
ylim([-0.1 0.5])
refline(0, nanmean(all.StabilityData.PDG.RDI_control(good_cells)));
refline(0, nanmean(all.StabilityData.NatMov.RDI_control(good_cells)));


% NOTE: statistical testing of RDI values on each session is handled by LMEM using the MixedEffectsModels_RDI function

%% responsiveness to stimuli proportions

% Add all fields to this object
all = StabilityAnalyzer;     
all.importData;         % import datafiles for first mouse ('RespData', 'RoiINFO', 'StabilityData')
all.addData();          % add datafiles for subsequent mice (all consecutive sessions)
all.addData([1 1 0 1 1 1]);         % if necessary, use input arg to indicate weeks for which data is present, only necessary for mice with non-consecutive sessions. 
                                    % e.g. [1 1 1 0 1 1 1] for data corresponding to sessions 1 2 3 5 6 7, [1 1 0 1 1 1] for sessions 1 2 4 5 6, etc.
all.setQualityThreshold(3);
all.cellSelection;                  % select 'quality' & 'presence'
qual = all.getUse_cells;

PDG_only = zeros(1, length(num_session_allmice));
NatMov_only = zeros(1, length(num_session_allmice));
both = zeros(1, length(num_session_allmice));
none = zeros(1, length(num_session_allmice));
% using pooled data
PDG_Responsive = all.RoiINFO.PDG_Responsive_thresh;  
NatMov_Responsive = all.RoiINFO.NatMov_Responsive_thresh;   

for ii = 1:length(num_session_allmice)
    try
        idx = cellstart_idx(ii):cellstart_idx(ii+1);
    catch
        idx = cellstart_idx(ii):length(qual);
    end
    PDG_only(ii) = sum(qual(idx) & PDG_Responsive(idx) & ~NatMov_Responsive(idx))/sum(qual(idx));
    NatMov_only(ii) = sum(qual(idx) & ~PDG_Responsive(idx) & NatMov_Responsive(idx))/sum(qual(idx));
    both(ii) = sum(qual(idx) & PDG_Responsive(idx) & NatMov_Responsive(idx))/sum(qual(idx));
    none(ii) = sum(qual(idx) & ~PDG_Responsive(idx) & ~NatMov_Responsive(idx))/sum(qual(idx));
end

mean_PDG_only = mean(PDG_only);
SEM_PDG_only = std(PDG_only, [], 2)/sqrt(length(PDG_only));
mean_NatMov_only = mean(NatMov_only);
SEM_NatMov_only = std(NatMov_only, [], 2)/sqrt(length(NatMov_only));
mean_both = mean(both);
SEM_both = std(both, [], 2)/sqrt(length(both));
mean_none = mean(none);
SEM_none = std(none, [], 2)/sqrt(length(none));
 
figure
pie([mean_PDG_only mean_NatMov_only mean_both mean_none]);


%% comparing response magnitudes between the two stimuli (dual responsive neurons)

nat_F0 = logical([zeros(1, 30) ones(1, 20) zeros(1, 300)]); 
pdg_F0 = logical(repmat([zeros(1, 20) ones(1, 20) zeros(1, 20)], 1, 12));

% uses same StabilityAnalyzer object as generated above
all.setQualityThreshold(3);
all.cellSelection;
num_sessions = all.num_sessions;

qual = all.getUse_cells;
sessions = all.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == num_sessions;
good_cells = qual' & present;

pdgResp = all.RespData.PDG.RespMat_Full(:, :, good_cells, :);
natResp = all.RespData.NatMov.RespMat_Full(:, :, good_cells, :);

pdg_catResp = [];
nat_catResp = [];
for kk = 1:all.num_sessions
    pdg_catResp = cat(1, pdg_catResp, pdgResp(:, :, :, kk));
    nat_catResp = cat(1, nat_catResp, natResp(:, :, :, kk));
end

mean_pdg_catResp = squeeze(nanmean(pdg_catResp, 1));
mean_nat_catResp = squeeze(nanmean(nat_catResp, 1));

z_pdg = mean(mean_pdg_catResp(~pdg_F0, :), 1)./std(mean_pdg_catResp(pdg_F0, :), [], 1);
z_nat = mean(mean_nat_catResp(~nat_F0, :), 1)./std(mean_nat_catResp(nat_F0, :), [], 1);

figure; axis square;
scatter(z_pdg, z_nat, 'filled');

[r, p] = corr(z_pdg', z_nat');



%% plotting response heatmaps and averages
vis = ResponseVisualizer();
vis.importData();               % import data for first field
vis.addData();                  % add data for additional fields (or just look at them one at a time)
vis.setQualityThreshold(3);     % set quality threshold
vis.cellSelection;              % select desired inclusion criteria
vis.plotPDGvNatMov_traces2(1, 'minimum');               % input args: starting neuron, normalization method ('minimum' is default)
vis.plotPDGvNatMov_avgs(23);            % input arg: starting neuron
vis.plotSingleStimulus_traces('NatMov', 1, 'first mode');

subplot(1, 2, 1)
set(gca, 'TickLength', [0 0]);
box off

%extract trial-averages for example neuron
PDG_mean = squeeze(mean(vis.RespData.PDG.RespMat_Full(:, :, 76, :), 1));
NatMov_mean = squeeze(mean(vis.RespData.NatMov.RespMat_onTime(:, :, 76, :), 1));


%% plotting polar tuning curves

%example neuron from a field
JCS003 = StabilityAnalyzer();
JCS003.importData();            % import 'RespData', 'RoiINFO', & 'OrientationData' datafiles
num_sessions = JCS003.num_sessions;
oriResp = squeeze(JCS003.OrientationData.oriResp(:, 76, :));
figure
theta = [0 2*pi];
pho = [20 20]; % set min
h1 = polarplot(theta,pho); % plot white point
set(h1,'color',[1 1 1])
theta = [0 2*pi];
hold on
num_ori = size(oriResp, 2);
upsample = 12;
sample1 = 1:upsample:num_ori*upsample;
sample2 = 1:1:num_ori*upsample;
theta = [0:360/size(oriResp, 2)/upsample:360]*pi/180;
pho = zeros(length(theta), num_sessions);
for kk = 1:num_sessions
    curr_pho = cat(2, interp1(cat(2, sample1, upsample*num_ori), cat(2, squeeze(oriResp(kk,:)), oriResp(kk, 1)), sample2, 'linear', 'extrap'), oriResp(kk,1))';
    curr_pho(curr_pho < 0) = 0; 
    h1 = polarplot(theta,curr_pho);
    pho(:, kk) = curr_pho;
end
% set(h1,'linewidth',2,'color',[0 0.5 0.8])
rlim([0 60])
title('Orientation tuning')


% pooled data
all = StabilityAnalyzer;
all.importData;                 % see above
all.addData;
all.addData([1 1 1 1 1 1 1]);

all.setQualityThreshold(3);
all.cellSelection;
good_cells = all.getUse_cells;
presence = all.getUse_sessions;
presence(isnan(presence)) = 0;
presence = logical(presence);
num_sessions = all.num_sessions;
num_cells = all.num_cells;

oriResp = all.OrientationData.oriResp;
oriPref = all.OrientationData.oriPref;
% circshift orientations based on prefOri
for ii = 1:num_cells
    currpref = oriPref(1, ii);
    shift = size(oriResp, 3) - currpref + 4;
    oriResp(:, ii, :) = circshift(squeeze(oriResp(:, ii, :)), shift, 2);
end

tuned = logical(all.OrientationData.isTuned(1, :));          % try with D0 tuned first
oriResp_tuned = oriResp(:, tuned & good_cells, :);
oriResp_avg = zeros(num_sessions, size(oriResp_tuned, 3));
for kk = 1:num_sessions
    oriResp_avg(kk, :) = nanmean(oriResp_tuned(kk, presence(tuned & good_cells, kk), :));
    n(kk) = sum(presence(tuned & good_cells, kk));
end

%percent tuned
isTuned = all.OrientationData.isTuned;
isTuned(isnan(isTuned)) = 0;
isTuned = logical(isTuned);
sessions = all.getUse_sessions;
for kk = 1:num_sessions
    fracTuned(kk) = sum(isTuned(kk, :) & good_cells & presence(:, kk)')/sum(good_cells);
end

% with upsampling
figure
theta = [0 2*pi];
pho = [20 20]; % set min
h1 = polarplot(theta,pho); % plot white point
set(h1,'color',[1 1 1])
theta = [0 2*pi];
hold on
num_ori = size(oriResp_avg, 2);
upsample = 12;
sample1 = 1:upsample:num_ori*upsample;
sample2 = 1:1:num_ori*upsample;
theta = [0:360/size(oriResp_avg, 2)/upsample:360]*pi/180;
pho = zeros(length(theta), num_sessions);
for kk = 1:num_sessions
    curr_pho = cat(2, interp1(cat(2, sample1, upsample*num_ori), cat(2, squeeze(oriResp_avg(kk,:)), oriResp_avg(kk, 1)), sample2, 'linear', 'extrap'), oriResp_avg(kk,1))';
    curr_pho(curr_pho<0) = 0; 
    h1 = polarplot(theta,curr_pho);
    pho(:, kk) = curr_pho;
end
% set(h1,'linewidth',2,'color',[0 0.5 0.8])
rlim([0 60])
title('Orientation tuning')

% % without upsampling
% theta = [0 2*pi];
% pho = [20 20]; % set min
% h1 = polarplot(theta,pho); % plot white point
% set(h1,'color',[1 1 1])
% theta = [0 2*pi];
% hold on
% num_ori = size(oriResp_avg, 2);
% theta = [0:360/size(oriResp_avg, 2):360]*pi/180;
% for kk = 1:num_sessions
%     pho = cat(2, squeeze(oriResp_avg(kk,:)), oriResp_avg(kk,1))';
%     pho(pho<0) = 0; 
%     h1 = polarplot(theta,pho);
% end
% % set(h1,'linewidth',2,'color',[0 0.5 0.8])
% title('Orientation tuning')




for ii = 2:all.num_cells

    theta = [0 2*pi];
    pho = [50 50]; % set min
    h1 = polarplot(theta,pho); % plot white point
    set(h1,'color',[1 1 1])
    theta = [0 2*pi];
    hold on
    num_ori = size(oriResp, 3);
    upsample = 12;
    sample1 = 1:upsample:num_ori*upsample;
    sample2 = 1:1:num_ori*upsample;
    theta = [0:360/size(oriResp, 3)/upsample:360]*pi/180;
    for kk = 1:all.num_sessions
        pho = cat(2, interp1(sample1, squeeze(oriResp(kk,ii,:)), sample2), squeeze(oriResp(kk,ii,1)))';
        pho(pho<0) = 0; % rectify
        h1 = polarplot(theta,pho);
    end
    set(h1,'linewidth',2,'color',[0 0.5 0.8])
    title('Orientation tuning')

    pause
    close

end









