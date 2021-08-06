function [A_cat, A_avg, trial_A_cat, gains_cat] = pupilsize_comparisons(week_vec, color, A_cat, A_avg, trial_A_cat, gains_cat)
% week_vec = week indicator for data, starting on session 1. e.g. [1 1 1 1 1 1] is sessions 1-6. [1 1 1 0 1 1 1] is sessions 1, 2, 3, 5, 6, 7
% color = field's color for scatter
%   A_cat = struct of previously built vectors of change in avg pupil size on each session for each stimulus (generated through previous iteration of function)
%   A_avg = struct of previously built vectors of change in avg pupil size over time for each stimulus (generated through previous iteration of function)
%   trial_A_cat = struct of previously built vectors of avg pupil size for each trial for each stimulus (generated through previous iteration of function)
%   gains_cat = struct of previously built vectors of by-trial DF/F gain factors for each stimulus (generated through previous iteration of function)


force_axes.pdg.query = 'No';               % change to Yes if you want to force data to certain axes (place data points outside limits on plot boundary)
force_axes.mov.query = 'Yes';
force_axes.x = [0 250];
force_axes.y = [0.4 1.6];

%% PDG vs MOV avg pupil size, each dot is one session
fprintf('Import eyetracking data.\n');
[pupil_fn, pupil_pn] = uigetfile('.mat');
cd(pupil_pn);
pupildata = importdata(pupil_fn);

num_sessions = length(pupildata.area);
pdg_A_vec = zeros(1, num_sessions);
mov_A_vec = zeros(1, num_sessions);

for kk = 1:num_sessions
    movarea = pupildata.area(kk).NatMov.raw;
    pdgarea = pupildata.area(kk).PDG.raw;
    
    if kk == 1
        ref_mov = mean(movarea);
        ref_pdg = mean(pdgarea);
    end

    mov_A_vec(kk) = (mean(movarea)/ref_mov)*100;                  % express as percentage of session 1 area
    pdg_A_vec(kk) = (mean(pdgarea)/ref_pdg)*100; 

end

% plot each mouse as different color

figure(1);
hold on
scatter(pdg_A_vec(2:end), mov_A_vec(2:end), 60, 'filled', 'MarkerFaceColor', color);

xlim([0 180])
ylim([0 180])
xlabel('PDG')
ylabel('MOV')
axis square
refline(1, 0)

% 2) pool data (for stat testing)
if nargin < 3
    A_cat.MOV = [];
    A_cat.PDG = [];
end

A_cat.MOV = [A_cat.MOV mov_A_vec(2:end)];
A_cat.PDG = [A_cat.PDG pdg_A_vec(2:end)];

% 3) test
[~, p, ~, stats] = ttest(A_cat.PDG', A_cat.MOV');
title(sprintf('PDG vs MOV avg. pupil area, p = %.2e', p));


%% Pupil area curves and RDI curves for each mouse

% Each mouse individually
% import data and produce vectors

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
plot(plot_weeks, pdg_A_vec, 'Color', 'r', 'LineStyle', '--');
plot(plot_weeks, mov_A_vec, 'Color', 'g', 'LineStyle', '--');
ylim([0 max([pdg_A_vec mov_A_vec])+50])
axis square
xlim([0 plot_weeks(end)+1])
hold off

%% average area curves across mice (all curves)


if nargin < 3
    A_avg.PDG = pdg_A_vec;
    A_avg.MOV = mov_A_vec;
else
    idx = 1;
    for kk = 1:length(week_vec)
        if week_vec(kk)
            try
                A_avg.PDG(kk) = mean([A_avg.PDG(kk) pdg_A_vec(idx)]);
                A_avg.MOV(kk) = mean([A_avg.MOV(kk) mov_A_vec(idx)]);
                idx = idx+1;
            catch
                A_avg.PDG = [A_avg.PDG pdg_A_vec(idx)];
                A_avg.MOV = [A_avg.MOV mov_A_vec(idx)];
            end
        end
    end
end

figure(4)
hold on
plot(plot_weeks, pdg_A_vec, 'Color', 'r', 'LineStyle', '--');          % plot individual mice in dotted lone
plot(plot_weeks, mov_A_vec, 'Color', 'g', 'LineStyle', '--');
ylim([0 max([pdg_A_vec mov_A_vec])+50])
axis square
xlim([0 8])
xlabel('Session')
ylabel('Eye movement (% of reference)')
title('Individual mice');

figure(5)
hold on
plot(1:length(A_avg.PDG), A_avg.PDG, 'Color', 'r', 'LineStyle', '-');                      % plot average as solid line
plot(1:length(A_avg.MOV), A_avg.MOV, 'Color', 'g', 'LineStyle', '-');
ylim([0 max([A_avg.PDG A_avg.MOV])+50])
axis square
xlim([0 8])
xlabel('Session')
ylabel('Eye movement (% of reference)')
hold off
title('Average across mice')


%% pupil size vs trial gain factor

A_bytrial_PDG = [];                                     % average pupil area on each trial expressed as a percentage of average pupil area across the session's recording
A_bytrial_MOV = [];

for kk = 1:num_sessions
    movarea = pupildata.area(kk).NatMov.sorted;
    pdgarea = pupildata.area(kk).PDG.sorted;
    
    A_bytrial_PDG = cat(1, A_bytrial_PDG, 100*mean(pdgarea, 2)/mean(mean(pdgarea, 1), 2));
    A_bytrial_MOV = cat(1, A_bytrial_MOV, 100*mean(movarea, 2)/mean(mean(movarea, 1), 2));
end

fprintf('Select threshold quality neurons responsive to MOV.\n');
an.cellSelection();             % select NatMov responsive for MOV
qual = an.getUse_cells;
sessions = an.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == num_sessions;
mov_cells = qual' & present;

fprintf('Select threshold quality neurons responsive to PDG.\n');
an.cellSelection();             % select PDG responsive for PDG
qual = an.getUse_cells;
sessions = an.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == num_sessions;
pdg_cells = qual' & present;

PDG_gains = zeros(1, 8*num_sessions);
MOV_gains = zeros(1, 30*num_sessions);
for kk = 1:num_sessions
    PDGresp = an.RespData.PDG.RespMat_Full(:, :, pdg_cells, kk);
    MOVresp = an.RespData.NatMov.RespMat_Full(:, :, mov_cells, kk);
    for tt = 1:size(PDGresp, 1)
        alphas = zeros(1, size(PDGresp, 3));
        for ii = 1:size(PDGresp, 3)
            curr_resp = PDGresp(:, :, ii);
            f = polyfit(squeeze(mean(curr_resp, 1)), curr_resp(tt, :), 1);
            alphas(ii) = f(1);
        end
        PDG_gains(8*(kk-1)+tt) = mean(alphas);              % either mean or median
    end   
    
    for tt = 1:size(MOVresp, 1)
        alphas = zeros(1, size(MOVresp, 3));
        for ii = 1:size(MOVresp, 3)
            curr_resp = MOVresp(:, :, ii);
            f = polyfit(squeeze(mean(curr_resp, 1)), curr_resp(tt, :), 1);
            alphas(ii) = f(1);
        end
        MOV_gains(30*(kk-1)+tt) = mean(alphas);              % either mean or median
    end
end



%pooled

if nargin < 3
    trial_A_cat.PDG = [];
    trial_A_cat.MOV = [];
    gains_cat.PDG = [];
    gains_cat.MOV = [];
end

% generate data for each mouse above, then add to vectors
trial_A_cat.MOV = cat(2, trial_A_cat.MOV, A_bytrial_MOV');
trial_A_cat.PDG = cat(2, trial_A_cat.PDG, A_bytrial_PDG');
gains_cat.MOV = [gains_cat.MOV MOV_gains];
gains_cat.PDG = [gains_cat.PDG PDG_gains];

[r_mov, p_mov] = corr(trial_A_cat.MOV', gains_cat.MOV');
[r_pdg, p_pdg] = corr(trial_A_cat.PDG', gains_cat.PDG');


% forcing axes limits for prettier plots
plot_trial_A_cat.PDG = trial_A_cat.PDG;
plot_gains_cat.PDG = gains_cat.PDG;
plot_trial_A_cat.MOV = trial_A_cat.MOV;
plot_gains_cat.MOV = gains_cat.MOV;

if strcmp(force_axes.pdg.query, 'Yes')
    plot_trial_A_cat.PDG(plot_trial_A_cat.PDG < force_axes.x(1)) = force_axes.x(1);
    plot_trial_A_cat.PDG(plot_trial_A_cat.PDG > force_axes.x(2)) = force_axes.x(2);
    plot_gains_cat.PDG(plot_gains_cat.PDG < force_axes.y(1)) = force_axes.y(1);
    plot_gains_cat.PDG(plot_gains_cat.PDG > force_axes.y(2)) = force_axes.y(2);  
end

if strcmp(force_axes.mov.query, 'Yes')
    plot_trial_A_cat.MOV(plot_trial_A_cat.MOV < force_axes.x(1)) = force_axes.x(1);
    plot_trial_A_cat.MOV(plot_trial_A_cat.MOV > force_axes.x(2)) = force_axes.x(2);
    plot_gains_cat.MOV(plot_gains_cat.MOV < force_axes.y(1)) = force_axes.y(1);
    plot_gains_cat.MOV(plot_gains_cat.MOV > force_axes.y(2)) = force_axes.y(2);
end

% plot
figure(6)
scatter(plot_trial_A_cat.MOV, plot_gains_cat.MOV, [], 'r', 'filled', 'MarkerFaceAlpha', 0.6);
hold on
xlabel('Pupil size (% of session-avg)')
ylabel('DFF trial gain factor')
title(sprintf('MOV, r = %.2f, p = %.2e', r_mov, p_mov))
axis square

figure(7)
scatter(plot_trial_A_cat.PDG, plot_gains_cat.PDG, [], 'b', 'filled', 'MarkerFaceAlpha', 0.6);
hold on
xlabel('Pupil size (% of session-avg)')
ylabel('DFF trial gain factor')
title(sprintf('PDG, r = %.2f, p = %.2e', r_pdg, p_pdg))
axis square

end