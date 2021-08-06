
%% First eyetracking supplementary figure: eye movements %%

function [WS_cat, WS_avg] = eyemovement_comparisons(week_vec, color, WS_cat, WS_avg)
% inputs:
%   week_vec = week indicator for field's data (e.g. [1 1 1 0 1 1 1] indicates data is from weeks 1, 2, 3, 5, 6, 7)
%   color = color for this field for scatter
%   WS_cat = struct of previously built vectors of change in avg distance from median centroid location on each session for each stimulus (generated through previous iteration of function)
%   WS_avg = struct of previously built vectors of change in avg eye movement over time for each stimulus (generated through previous iteration of function)
%
% outputs:
%   WS_cat = see above, pooled data structure built with successive calls of the function using the input arg
%   WS_avg = see above, pooled data structure built with successive calls of the function using the input arg
%
% see master script for generation of figures


%% PDG vs MOV WS movement scatterplot
% avg distance from median centroid location on each session, PDG vs MOV
% do this separately to create data vectors for each stimulus
% 1) extract data and create data vector for each mouse
fprintf('Import eyetracking data.\n');
[pupil_fn, pupil_pn] = uigetfile('.mat');
cd(pupil_pn);
pupildata = importdata(pupil_fn);

num_sessions = length(pupildata.centroid);
pdg_WS_vec = zeros(1, num_sessions);
mov_WS_vec = zeros(1, num_sessions);

for kk = 1:num_sessions
    movcentroid = pupildata.centroid(kk).NatMov.raw;
    pdgcentroid = pupildata.centroid(kk).PDG.raw;
    movmed = median(movcentroid);
    pdgmed = median(pdgcentroid);
    
    movdis = zeros(1, size(movcentroid, 1));
    for tt = 1:size(movcentroid, 1)
        movdis(tt) = pdist([movcentroid(tt, :); movmed], 'euclidean');
    end
    curr_dis_mov = mean(movdis);
    
    pdgdis = zeros(1, size(pdgcentroid, 1));
    for tt = 1:size(pdgcentroid, 1)
        pdgdis(tt) = pdist([pdgcentroid(tt, :); pdgmed], 'euclidean');
    end
    curr_dis_pdg = mean(pdgdis);
    
    if kk == 1
        ref_mov = curr_dis_mov;
        ref_pdg = curr_dis_pdg;
    end
    
    mov_WS_vec(kk) = (mean(movdis)/ref_mov)*100;                     % this is expressed as percentage of reference session value
    pdg_WS_vec(kk) = (mean(pdgdis)/ref_pdg)*100;

end

% plot each mouse as different color

figure(1);    
scatter(pdg_WS_vec(2:end), mov_WS_vec(2:end), 60, 'filled', 'MarkerFaceColor', color);     
hold on


xlabel('PDG')
ylabel('MOV')
xlim([0 200])
ylim([0 200])
axis square
% refline(1, 0)

% 2) pool data (for stat testing)
% initialize vectors once for the first mouse input
if nargin < 3
    WS_cat.MOV = [];
    WS_cat.PDG = [];
end

% run for each mouse
WS_cat.MOV = [WS_cat.MOV mov_WS_vec(2:end)];
WS_cat.PDG = [WS_cat.PDG pdg_WS_vec(2:end)];

% 3) test

% [r, p] = corr(WS_cat.PDG', WS_cat.MOV');
% refline(1, 0)

[~, p, ~, stats] = ttest(WS_cat.PDG', WS_cat.MOV');
title(sprintf('Avg distance from med. centroid loc., p = %.2e', p));

%% WS eye movement curves and RDI curves for each mouse
plot_weeks = find(week_vec);

an = StabilityAnalyzer();
fprintf('Load in datafiles (RespData, RoiINFO, StabilityData)\n');
an.importData();
an.setQualityThreshold(3);
figure(2)
fprintf('Plotting RDIcurves to extract averages. Select inclusion criteria\n');
an.RDIplotter('PDG', 'Average');
hold on
an.RDIplotter('NatMov', 'Average');
RDIdata = an.getRDIdata;
pdgRDI = RDIdata.PDG.RDI_avg;
movRDI = RDIdata.NatMov.RDI_avg;
pdgRDIsem = RDIdata.PDG.RDI_SEM;
movRDIsem = RDIdata.NatMov.RDI_SEM;

figure(3)
hold on
yyaxis left
errorbar(plot_weeks(2:end), pdgRDI, pdgRDIsem, 'Color', 'r', 'LineStyle', '-');
errorbar(plot_weeks(2:end), movRDI, movRDIsem, 'Color', 'g', 'LineStyle', '-');
ylim([0 0.45])
yyaxis right
plot(plot_weeks, pdg_WS_vec, 'Color', 'r', 'LineStyle', '--');
plot(plot_weeks, mov_WS_vec, 'Color', 'g', 'LineStyle', '--');
ylim([0 max([pdg_WS_vec mov_WS_vec])+50])
axis square
xlim([0 plot_weeks(end)+1])
hold off

%% average WS mov curves across mice (all curves)

% start with SKKS066
if nargin < 3
    WS_avg.PDG = pdg_WS_vec;
    WS_avg.MOV = mov_WS_vec;
else
    idx = 1;
    for kk = 1:length(week_vec)
        if week_vec(kk)
            try
                WS_avg.PDG(kk) = mean([WS_avg.PDG(kk) pdg_WS_vec(idx)]);
                WS_avg.MOV(kk) = mean([WS_avg.MOV(kk) mov_WS_vec(idx)]);
                idx = idx+1;
            catch
                WS_avg.PDG = [WS_avg.PDG pdg_WS_vec(idx)];
                WS_avg.MOV = [WS_avg.MOV mov_WS_vec(idx)];
            end
        end
    end
end
        

figure(4)
hold on
% yyaxis left
% errorbar(2:7, pdgRDI, pdgRDIsem, 'Color', 'r', 'LineStyle', '-');
% errorbar(2:7, movRDI, movRDIsem, 'Color', 'g', 'LineStyle', '-');
% ylim([0 0.4])
% yyaxis right
plot(plot_weeks, pdg_WS_vec, 'Color', 'r', 'LineStyle', '--');          % plot individual mice in dotted lone
plot(plot_weeks, mov_WS_vec, 'Color', 'g', 'LineStyle', '--');
ylim([0 max([pdg_WS_vec mov_WS_vec])+50])
axis square
xlim([0 8])
xlabel('Session')
ylabel('Eye movement (% of reference)')
title('Individual mice');

figure(5)
hold on
plot(1:length(WS_avg.PDG), WS_avg.PDG, 'Color', 'r', 'LineStyle', '-');                      % plot average as solid line
plot(1:length(WS_avg.MOV), WS_avg.MOV, 'Color', 'g', 'LineStyle', '-');
ylim([0 max([WS_avg.PDG WS_avg.MOV])+50])
axis square
xlim([0 8])
xlabel('Session')
ylabel('Eye movement (% of reference)')
hold off
title('Average across mice')
end



