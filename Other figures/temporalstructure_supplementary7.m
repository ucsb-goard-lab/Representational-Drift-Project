% temporal structure figurescript

% For RDI curves: use LME function. This is just for generating the control lines for the plot
temp = StabilityAnalyzer();
temp.importData();                      % import first field
temp.addData();                         % add subsequent fields (use input arg to indicate weeks if necessary)
temp.setQualityThreshold(3);
temp.cellSelection;
qual = temp.getUse_cells;
sessions = temp.getUse_sessions;
sessions(isnan(sessions)) = 1;            % nan values converted to zeros for logical purposes
present = sum(sessions, 2) == temp.num_sessions;
good_cells = qual' & present;




a = refline(0, nanmean(temp.StabilityData.PDG.RDI_control(good_cells)));    a.Color = 'b';
b = refline(0, nanmean(temp.StabilityData.NatMov.RDI_control(good_cells))); b.Color = 'r';
c = refline(0, nanmean(temp.StabilityData.PDGcont.RDI_control(good_cells)));    c.Color = 'g';
d = refline(0, nanmean(temp.StabilityData.NatMovdisc.RDI_control(good_cells))); d.Color = 'k';