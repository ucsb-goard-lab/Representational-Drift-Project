% inhibitory populations


%% Example neuron images
gad = ChronicImaging_Mouse;
gad.importData();           % import data for example field
gad.evaluateROIs;           % run 'Visual_Inspection_A' to process cell images

num_sessions = gad.num_sessions;
neuron = 81;
map_type = 'NatMov';
% cell maps
figure
for ii = neuron:gad.num_cells % ii = neuron;
    for kk = 1:num_sessions
        subplot(2, 7, kk)
        imagesc(imgaussfilt(gad.Vis_cellInspection_maps.(map_type).cell_avgproj_mat(:, :, kk, ii)));
        colormap gray
        axis square
    end
    for kk = 1:num_sessions
        subplot(2, 7, kk+7)
        imagesc(imgaussfilt(gad.Vis_cellInspection_maps.(map_type).cell_actmap_mat(:, :, kk, ii)));
        colormap gray
        axis square
    end
    title(sprintf('neuron %d', ii))
    pause
end

%% field image
load('TSeries-09202019-1437-NatMov-002_registered_data.mat')
imagesc(data.avg_projection); colormap gray; axis square


avgproj = data.avg_projection;
MIN = min(min(avgproj));
avgproj = avgproj - MIN;
MAX = max(max(avgproj));
avgproj = avgproj/MAX;

adjusted = imadjust(avgproj);
imagesc(adjusted); colormap gray; axis square;
set(gca, 'xtick', []);
set(gca, 'ytick', []);


%% RDI curves
gad = StabilityAnalyzer();
% import first field
gad.importData();   
% add subsequent fields
gad.addData();                  % use input arg to indicate weeks if a mouse has data on non-consecutive weeks (e.g. [1 1 1 0 1 1 1] indicates data corresponds to weeks 1, 2, 3, 5, 6, 7)
gad.setQualityThreshold(3);
gad.cellSelection();
good_cells = gad.getUse_cells;

figure
gad.RDIplotter('PDG', 'Average');
hold on
gad.RDIplotter('NatMov', 'Average');
ylim([-0.1 0.5]);
legend off
refline(0, nanmean(gad.StabilityData.NatMov.RDI_control(good_cells)));
refline(0, nanmean(gad.StabilityData.PDG.RDI_control(good_cells)));

% NOTE: RDI curves shown on figure and stats are handled by LMEM using the MixedEffectsModels_RDI function


%responsivity
gad.setQualityThreshold(3);         % uses StabilityAnalyzer object with pooled data generated above
gad.cellSelection;      % select quality + presence
qual = gad.getUse_cells;
sessions = gad.getUse_sessions;
sessions(isnan(sessions)) = 0;            % nan values converted to zeros for logical purposes
good_cells = qual' & sessions;                % cells of sufficient quality & presence on any given session
D0_good_cells = good_cells(:, 1);                             % quality cells present on D0

gad.cellSelection;               % select quality + all responsiveness criteria
dual_responsive = gad.getUse_cells';       %cells responsive to both stimuli

num_session_allmice = [7 7 7 7];                %number of sessions for each mouse in the order in which they were loaded in to the Analyzer above
num_cells_allmice = gad.getCellList;

ct = 1;
cellstart_idx = [];
for ii = 1:length(num_session_allmice)
    cellstart_idx(ii) = ct;                                    % starting index for the first cell in every mouse
    ct = ct + num_cells_allmice(ii);
end

PDG_only = zeros(1, length(num_session_allmice));
NatMov_only = zeros(1, length(num_session_allmice));
both = zeros(1, length(num_session_allmice));
none = zeros(1, length(num_session_allmice));

PDG_Responsive = D0_good_cells' & gad.RoiINFO.PDG_Responsive_thresh; 
NatMov_Responsive = D0_good_cells' & gad.RoiINFO.NatMov_Responsive_thresh;    

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


% single neuron RDI curve
neuron = 229;
NatMovcurve = gad.StabilityData.NatMov.RDI(:, neuron);
PDGcurve = gad.StabilityData.PDG.RDI(:, neuron);
figure
axis square
hold on
plot(PDGcurve, 'color', [0 0.7 0.7])
plot(NatMovcurve, 'color', [0.8 0.5 0.5])
ylim([-0.1 0.5])
xlim([0 gad.num_sessions])
refline(0, gad.StabilityData.NatMov.RDI_control(neuron))
refline(0, gad.StabilityData.PDG.RDI_control(neuron))


% neuron responses

gadvis = ResponseVisualizer();
gadvis.importData();
gadvis.addData();
gadvis.setQualityThreshold(3);
gadvis.cellSelection();
good_cells = gadvis.getUse_cells;

gadvis.plotPDGvNatMov_traces2(1, 'minimum')
box off
set(gca, 'xtick', [])
set(gca, 'ytick', [])

gadvis.plotPDGvNatMov_avgs(229, [0 0.45]);

%extract trial-averages for example neuron
PDG_mean = squeeze(mean(gadvis.RespData.PDG.RespMat_Full(:, :, 229, :), 1));
NatMov_mean = squeeze(mean(gadvis.RespData.NatMov.RespMat_onTime(:, :, 229, :), 1));

