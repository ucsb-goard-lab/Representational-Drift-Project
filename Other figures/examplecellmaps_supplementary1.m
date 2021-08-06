% initialize Analyzer
all = StabilityAnalyzer();
% First mouse: import data (RespData, RoiINFO, StabilityData)
all.importData();
% For subsequent mice, add data
all.addData();
% if the mouse has non-consecutive sessions, indicate weeks for which the mouse has data in the following way
all.addData([1 1 0 1 1 1]);             % indicates data present on weeks 1, 2, 4, 5, 6


% distribution of quality values
qual = all.RoiINFO.quality;
figure
histogram(qual);
xline(2.5)      % threshold indicator (2.5 == 3)
xlim([0 6])
axis square


%% Generate RDI curves using different quality thresholds

% Note: all pooled-data RDI curves are generated using the MixedEffectsModels_RDI function

% Generate control lines for RDI plots
all.setQualityThreshold(3);
all.cellSelection;
good_cells = all.getUse_cells;
refline(0, nanmean(all.StabilityData.PDG.RDI_control(good_cells)));
refline(0, nanmean(all.StabilityData.NatMov.RDI_control(good_cells)));


