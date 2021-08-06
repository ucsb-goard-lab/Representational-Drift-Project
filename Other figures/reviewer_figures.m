% REVIEWER FIGURES


%% age at first recording vs. final RDI
% set age, create analysis and import data for each mouse
age = 250;
an = StabilityAnalyzer();
an.importData();
% an.addData();                 % for mice with mulitple fields, combine into same object by adding subsequent fields
an.setQualityThreshold(3);


figure
an.RDIplotter('PDG', 'Average');         % running RDIplotter produces the RDIdata structure, which we can then extract
hold on
an.RDIplotter('NatMov', 'Average');

RDIdata = an.getRDIdata;
PDG_RDIavg = RDIdata.PDG.RDI_avg(end);
NatMov_RDIavg = RDIdata.NatMov.RDI_avg(end);

% initialize once
ages = [];
PDG_RDIs = [];
MOV_RDIs = [];

% find values for each mouse (above), then concatenate
ages = [ages age];
PDG_RDIs = [PDG_RDIs PDG_RDIavg];
MOV_RDIs = [MOV_RDIs NatMov_RDIavg];

% plot pooled data
figure; 
scatter(ages, PDG_RDIs, 60, 'filled');
hold on
scatter(ages, MOV_RDIs, 60, 'filled');
xlabel('Age (days)');
ylabel('Final RDI value')
[r_pdg, p_pdg] = corr(ages', PDG_RDIs');
[r_mov, p_mov] = corr(ages', MOV_RDIs'); 
title(sprintf('PDG r = %.2f, p = %.2f; MOV r = %.2f, p = %.2f', r_pdg, p_pdg, r_mov, p_mov));
axis square



















