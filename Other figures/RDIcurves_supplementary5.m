
an = StabilityAnalyzer;
an.importData();
% an.addData();
an.importData([1 1 1 1 1 1 1]);        % indicates weeks to which data corresponds. e.g. [1 1 1 0 1 1 1] indicates data corresponds to weeks 1, 2, 3, 5, 6, 7. If data is on all consecutive weeks, don't need this input argument
an.setQualityThreshold(3);
an.cellSelection;                      % select dual responsive neurons
weeks = [1 1 1 0 1 1 1];
an.setWeeks(weeks);

qual = an.getUse_cells;
sessions = an.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == an.num_sessions;
good_cells = qual' & present;


figure;
an.RDIplotter('PDG', 'Average');
hold on
an.RDIplotter('NatMov', 'Average');
legend off
ylim([-0.2 1])
a = refline(0, nanmean(an.StabilityData.PDG.RDI_control(good_cells)));  a.Color = 'b';
b = refline(0, nanmean(an.StabilityData.NatMov.RDI_control(good_cells))); b.Color = 'r';


RDIdata = an.getRDIdata;
PDG_RDI = RDIdata.PDG.RDI_included_scatter;
NatMov_RDI = RDIdata.NatMov.RDI_included_scatter;

p = zeros(1, an.num_sessions-1);
for kk = 1:size(PDG_RDI, 1)
    curr_pdg_vals = PDG_RDI(kk, :);
    curr_pdg_vals(isnan(curr_pdg_vals)) = [];
    curr_natmov_vals = NatMov_RDI(kk, :);
    curr_natmov_vals(isnan(curr_natmov_vals)) = [];
    [~,p(kk)] = ttest(curr_pdg_vals, curr_natmov_vals, 'tail', 'both');
end


titlestring = [];
for kk = 1:length(p)
    titlestring = sprintf('%s %0.3e', titlestring, p(kk));
end

title(titlestring);


% Adata.n11.PDG.curve = RDIdata.PDG.RDI_avg;
% Adata.n11.PDG.errorbars = RDIdata.PDG.RDI_SEM;
% Adata.n11.NatMov.curve = RDIdata.NatMov.RDI_avg;
% Adata.n11.NatMov.errorbars = RDIdata.NatMov.RDI_SEM;
% Adata.n11.PDG.control = nanmean(STDM077.StabilityData.PDG.RDI_control(good_cells));
% Adata.n11.NatMov.control = nanmean(STDM077.StabilityData.NatMov.RDI_control(good_cells));


%% slope and y-intercept
an = StabilityAnalyzer();
an.importData([1 1 1 0 1 1 1]);        % see above
% STDM077.addData();
an.setQualityThreshold(3);
an.cellSelection;
qual = an.getUse_cells;
num_cells = an.num_cells;

f = figure;
an.RDIplotter('PDG', 'Average');
hold on
an.RDIplotter('NatMov', 'Average');
legend off
ylim([-0.2 1])
RDIdata = an.getRDIdata();

week_vector = [1 2 3 4 5 6 7];      % [1:n] where n = last week for which mouse was imaged, regardless of non-consecutive sessions. e.g. mouse with data on session 1, 2, 4, 5, 6 is still 1, 2, 3, 4, 5, 6.
PDGavg = RDIdata.PDG.RDI_avg;
NatMovavg = RDIdata.NatMov.RDI_avg;
PDGall = RDIdata.PDG.RDI_included_scatter;
NatMovall = RDIdata.NatMov.RDI_included_scatter;

%for individual fields (comparing single cells)
pdg_ms = [];
nat_ms = [];
for ii = 1:num_cells
    if qual(ii)
        nanidx = isnan(PDGall(:, ii));
        currfit = polyfit(week_vector(~nanidx), PDGall(~nanidx, ii)', 1);
        pdg_ms = [pdg_ms currfit(1)];
        nanidx = isnan(NatMovall(:, ii));
        currfit = polyfit(week_vector(~nanidx), NatMovall(~nanidx, ii)', 1);
        nat_ms = [nat_ms currfit(1)];
    end
end

% remove outliers (only for plotting)
pdg_outliers = isoutlier(pdg_ms, 'mean');
nat_outliers = isoutlier(nat_ms, 'mean');

figure
hold on
for jj = 1:length(pdg_ms)
    if ~pdg_outliers(jj) && ~nat_outliers(jj)
        plot(1:2, [pdg_ms(jj) nat_ms(jj)], 'k');
        scatter(1, pdg_ms(jj), 'filled', 'MarkerFaceColor', 'b');
        scatter(2, nat_ms(jj), 'filled', 'MarkerFaceColor', 'r');
    end
end
plot(1:2, [mean(pdg_ms) mean(nat_ms)]', 'r');
scatter(1, mean(pdg_ms), 'filled', 'MarkerFaceColor', 'k');
scatter(2, mean(nat_ms), 'filled', 'MarkerFaceColor', 'k');

xlim([0 3])
ylim([-0.2 0.5])
[~, p, ~, stats] = ttest(pdg_ms, nat_ms);
title(sprintf('slope comparison, p = %.2e, n = %.2d', p, length(pdg_ms)));
set(gcf, 'Position', [400 400 300 800])

Adata.n8.PDG.slopes = pdg_ms;
Adata.n8.PDG.meanslope = mean(pdg_ms);
Adata.n8.NatMov.slopes = nat_ms;
Adata.n8.NatMov.meanslope = mean(nat_ms);

% for pooling
pdgfit = polyfit(week_vector, PDGavg, 1);
natfit = polyfit(week_vector, NatMovavg, 1);

% initialize once
pdg_m_cat = [];
natmov_m_cat = [];
pdg_b_cat = [];
natmov_b_cat = [];

% calculate values for each mouse (above), and then concatenate
pdg_m_cat = [pdg_m_cat pdgfit(1)];
natmov_m_cat = [natmov_m_cat natfit(1)];
pdg_b_cat = [pdg_b_cat pdgfit(2)];
natmov_b_cat = [natmov_b_cat natfit(2)];

save pdg_m_cat pdg_m_cat
save pdg_b_cat pdg_b_cat
save natmov_m_cat natmov_m_cat
save natmov_b_cat natmov_b_cat

% slopes
figure
hold on
for jj = 1:length(pdg_m_cat)
    plot(1:2, [pdg_m_cat(jj) natmov_m_cat(jj)]);
    scatter(1, pdg_m_cat(jj), 'filled');
    scatter(2, natmov_m_cat(jj), 'filled');
end
xlim([0 3])
ylim([-0.01 0.2])
[~, p, ~, stats] = ttest(pdg_m_cat, natmov_m_cat);
tstat = stats.tstat;
title(sprintf('slope comparison, p = %.2e', p));

% y-intercepts
figure
hold on
for jj = 1:length(pdg_b_cat)
    plot(1:2, [pdg_b_cat(jj) natmov_b_cat(jj)]);
    scatter(1, pdg_b_cat(jj), 'filled');
    scatter(2, natmov_b_cat(jj), 'filled');
end
xlim([0 3])
ylim([-0.2 0.25])
[~, p, ~, stats] = ttest(pdg_b_cat, natmov_b_cat);
tstat = stats.tstat;
title(sprintf('y-int comparison, p = %.2e', p));


