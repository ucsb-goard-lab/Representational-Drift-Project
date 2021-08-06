

an = StabilityAnalyzer();           % initialize object
an.importData();                    % import first field
an.addData();                       % add subsequent fields
an.setQualityThreshold(3);
an.cellSelection();                 % select 'quality', 'presence', 'PDG_Responsive_thresh', 'NatMov_Responsive_thresh'
qual = an.getUse_cells;
sessions = an.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == an.num_sessions;
good_cells = qual' & present;

% NOTE: RDI curves and stats in figure is handled by LMEM using MixedEffectsModels_RDI function
figure
an.RDIplotter('PDG', 'Average');
hold on
an.RDIplotter('MultiMov', 'Average');
ylim([-0.1 0.5])
legend off
% for control lines
a = refline(0, nanmean(an.StabilityData.PDG.RDI_control(good_cells)));    a.Color = 'b';
b = refline(0, nanmean(an.StabilityData.MultiMov.RDI_control{1}(good_cells)));   b.Color = 'r';
c = refline(0, nanmean(an.StabilityData.MultiMov.RDI_control{2}(good_cells)));    c.Color = [1 0.5 0];
d = refline(0, nanmean(an.StabilityData.MultiMov.RDI_control{3}(good_cells)));   d.Color = 'y';

