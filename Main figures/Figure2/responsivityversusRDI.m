all = StabilityAnalyzer;
all.importData();               % import first field
all.addData();                  % add subsequent fields
all.setQualityThreshold(3);
all.cellSelection();


corr_type = 'Pearson';
stimtype = 'NatMov';
timespan = 'average';

binsize = 10;        % percent


%% responsivity vs RDI
close all
[plotdata, stats] = all.responsivity('responsivity', stimtype, timespan, binsize);   

ylim([-0.12 0.8])
xlim([0 1])



