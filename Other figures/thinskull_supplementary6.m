thinskull = StabilityAnalyzer();
thinskull.importData();                 % import data for first field
thinskull.addData();                    % add data for subsequent fields
thinskull.setQualityThreshold(3);
thinskull.cellSelection;
num_sessions = thinskull.num_sessions;

qual = thinskull.getUse_cells;
sessions = thinskull.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == num_sessions;
good_cells = qual' & present;

% NOTE: RDI curves and stats are handled by LMEM in the MixedEffectsModels_RDI function
% figure
% thinskull.RDIplotter('PDG', 'Average');
% hold on
% thinskull.RDIplotter('NatMov', 'Average');
% ylim([-0.1 0.5])

%Generate control lines
refline(0, nanmean(thinskull.StabilityData.PDG.RDI_control(good_cells)));
refline(0, nanmean(thinskull.StabilityData.NatMov.RDI_control(good_cells)));

