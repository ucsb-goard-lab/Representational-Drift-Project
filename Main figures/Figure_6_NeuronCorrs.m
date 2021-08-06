
% create Analyzer (data import is in constructor in this one)
MJG013w1 = EnsembleAnalyzer();
% set quality threshold and assign inclusion criteria
MJG013w1.setQualityThreshold(3);
MJG013w1.cellSelection;
num_sessions = MJG013w1.num_sessions;
num_cells = MJG013w1.num_cells;

qual = MJG013w1.getUse_cells;            % reflects included neurons based on above selections
sessions = MJG013w1.getUse_sessions;     % presence data for each neuron
% consider neurons present on all sessions (some have NaN values where data was not collected, set to 1 for logical purposes)
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == num_sessions;
good_cells = qual' & present;           % included neurons that are present on all sessions

% if using same neurons for both stimuli, these will be identical
% if not, run lines 5-15 separately for each stim, and save the good_cells vector separately for each stimulus
save good_cells_natmov good_cells
save good_cells_pdg good_cells


%% ------- Signal corrs -----------%

% Full matrices
Scorr_natmov = zeros(num_cells, num_cells, num_sessions);
Scorr_pdg = zeros(num_cells, num_cells, num_sessions);
for kk = 1:num_sessions
    fprintf('Session %d of %d\n', kk, num_sessions);
    Scorr_natmov(:, :, kk) = MJG013w1.SignalCorr(MJG013w1.RespData.NatMov.RespMat_Full(:, :, :, kk));
    Scorr_pdg(:, :, kk) = MJG013w1.SignalCorr(MJG013w1.RespData.PDG.RespMat_Full(:, :, :, kk));
end
save Scorr_natmov Scorr_natmov
save Scorr_pdg Scorr_pdg


%% ------- Noise corrs 201210 ---------%

% Full matrices
Ncorr_natmov = zeros(sum(good_cells), sum(good_cells), num_sessions);           % need to limit to only good_cells to save time
Ncorr_pdg = zeros(sum(good_cells), sum(good_cells), num_sessions);



frames = 5760;      % to match amount of data used for corr calculation between stimuli
for kk = 1:num_sessions
    fprintf('Session %d out of %d\n', kk, num_sessions);
    PDG = MJG013w1.RespData.PDG.RespMat_Full(:, :, :, kk);               
    NatMov = MJG013w1.RespData.NatMov.RespMat_Full(:, :, :, kk);
    
    fprintf('PDG...\n');
    Ncorr_pdg(:, :, kk) = MJG013w1.NoiseCorr(PDG, frames);       % the NoiseCorr method only performs calculations on good_cells, as established above.
    
    fprintf('NatMov...\n');
    Ncorr_natmov(:, :, kk) = MJG013w1.NoiseCorr(NatMov, frames);

end

save Ncorr_natmov Ncorr_natmov
save Ncorr_pdg Ncorr_pdg

%% Visualization

% load either Scorrs or Ncorrs
corr_natmov = importdata('Ncorr_natmov.mat');
corr_pdg = importdata('Ncorr_pdg.mat');

% rectify identity line to 0 for visualization
num_cells = size(corr_natmov, 1);
num_sessions = size(corr_natmov, 3);
I = logical(eye(num_cells));
for kk = 1:num_sessions                         
    natmov_curr = corr_natmov(:, :, kk);
    natmov_curr(I) = 0;
    corr_natmov(:, :, kk) = natmov_curr;
    
    pdg_curr = corr_pdg(:, :, kk);
    pdg_curr(I) = 0;
    corr_pdg(:, :, kk) = pdg_curr;
end

% extract only good_cells -- only necessary for Scorrs, since Ncorrs are output with good_cells only already
load good_cells_natmov
corr_natmov = corr_natmov(good_cells, good_cells, :);
load good_cells_pdg
corr_pdg = corr_pdg(good_cells, good_cells, :);

% sort
PDGResp = MJG013w1.RespData.PDG.RespMat_Full(:, :, good_cells, :);
NatResp = MJG013w1.RespData.NatMov.RespMat_Full(:, :, good_cells, :);
PDGcatResp = [];
NatcatResp = [];
for kk = 1:MJG013w1.num_sessions
    PDGcatResp = cat(1, PDGcatResp, PDGResp(:, :, :, kk));
    NatcatResp = cat(1, NatcatResp, NatResp(:, :, :, kk));
end
[~, pdgmaxidcs] = max(squeeze(mean(PDGcatResp, 1)), [], 1);
[~, natmaxidcs] = max(squeeze(mean(NatcatResp, 1)), [], 1);

[~, pdgorder] = sort(pdgmaxidcs);
[~, natorder] = sort(natmaxidcs);

corr_pdg_sorted = corr_pdg(pdgorder, pdgorder, :);
corr_natmov_sorted = corr_natmov(natorder, natorder, :);


%%plotting

% A: plot corr matrices D0 and D42 for both stimuli
coloraxis = ([-0.2 0.2]);
figure
colormap(redblue);
subplot(2, 2, 1)
imagesc(corr_pdg_sorted(:, :, 1))
axis square
caxis(coloraxis)
colorbar
xlabel('D0')
ylabel('D0')
subplot(2, 2, 2)
imagesc(corr_pdg_sorted(:, :, num_sessions))
axis square
caxis(coloraxis)
colorbar
xlabel('D final')
ylabel('D final')
subplot(2, 2, 3)
imagesc(corr_natmov_sorted(:, :, 1))
axis square
caxis(coloraxis)
colorbar
xlabel('D0')
ylabel('D0')
subplot(2, 2, 4)
imagesc(corr_natmov_sorted(:, :, num_sessions))
axis square
caxis(coloraxis)
colorbar
xlabel('D final')
ylabel('D final')

% difference matrix between D0 and D42 for both stimuli
coloraxis = ([-1 1]);
figure
colormap(redblue);
subplot(1, 2, 1)
imagesc(corr_pdg_sorted(:, :, num_sessions) - corr_pdg_sorted(:, :, 1));
axis square
caxis(coloraxis)
colorbar
subplot(1, 2, 2)
imagesc(corr_natmov_sorted(:, :, num_sessions) - corr_natmov_sorted(:, :, 1));
axis square
caxis(coloraxis)
colorbar


%% changes in connectivity magnitude

%Mean connectivity
pdgdiff = zeros(1, sum(good_cells));
natdiff = zeros(1, sum(good_cells));
for ii = 1:sum(good_cells)
    pdgdiff(ii) = mean(corr_pdg(1:end~=ii, ii, num_sessions), 1) - mean(corr_pdg(1:end~=ii, ii, 1), 1);         % D42 mean conn. - D0 mean conn.
    natdiff(ii) = mean(corr_natmov(1:end~=ii, ii, num_sessions), 1) - mean(corr_natmov(1:end~=ii, ii, 1), 1);
end

% % B: Example field - distribution of neurons' D42 mean conn. - D0 mean conn., PDG vs NatMov
figure
histogram(abs(pdgdiff), 0:0.005:0.08);            % bins here are for Scorrs
hold on
histogram(abs(natdiff), 0:0.005:0.08);
xline(mean(abs(pdgdiff)));
xline(mean(abs(natdiff)));
xlabel('|Delta mean conn.|')
ylabel('Neurons')
[p, ~, stats] = ranksum(abs(pdgdiff), abs(natdiff));
title(sprintf('|D42 mean conn. - D0 mean conn.|, PDG v NatMov, p = %.2e', p))

% % C: pooling data for all fields - distribution of fields' mean [D42 mean conn. - D0 mean conn.] for PDG vs NatMov
% initialize these vectors once
pdg_meanconn_absdiff_cat = [];
nat_meanconn_absdiff_cat = [];

% % add each mouse 
pdg_meanconn_absdiff_cat = [pdg_meanconn_absdiff_cat mean(abs(pdgdiff))];          %mean of the distribution of connectivity difference abs values for each of the stimuli
nat_meanconn_absdiff_cat = [nat_meanconn_absdiff_cat mean(abs(natdiff))];

save pdg_meanconn_absdiff_cat pdg_meanconn_absdiff_cat
save nat_meanconn_absdiff_cat nat_meanconn_absdiff_cat


%% Correlating corr matrices

% CCbs

%import either Scorr or Ncorr
corr_natmov = importdata('Ncorr_natmov.mat');
corr_pdg = importdata('Ncorr_pdg.mat');

%rectify identity line to 0 
num_cells = size(corr_natmov, 1);
num_sessions = size(corr_natmov, 3);
I = logical(eye(num_cells));
for kk = 1:num_sessions                         
    natmov_curr = corr_natmov(:, :, kk);
    natmov_curr(I) = 0;
    corr_natmov(:, :, kk) = natmov_curr;
    
    pdg_curr = corr_pdg(:, :, kk);
    pdg_curr(I) = 0;
    corr_pdg(:, :, kk) = pdg_curr;
end

%extract only good_cells if necessary
load good_cells_natmov
corr_natmov = corr_natmov(good_cells, good_cells, :);
load good_cells_pdg
corr_pdg = corr_pdg(good_cells, good_cells, :);

nat_btwses_cc = zeros(1, num_sessions-1);
pdg_btwses_cc = zeros(1, num_sessions-1);

%2D correlation
for kk = 2:num_sessions
    nat_btwses_cc(kk-1) = corr2(corr_natmov(:, :, kk), corr_natmov(:, :, 1));
    pdg_btwses_cc(kk-1) = corr2(corr_pdg(:, :, kk), corr_pdg(:, :, 1));
end
% gives the same result
% for kk = 2:num_sessions
%     nat_btwses_cc(kk-1) = max(max(normxcorr2(corr_natmov(:, :, kk), corr_natmov(:, :, 1))));
%     pdg_btwses_cc(kk-1) = max(max(normxcorr2(corr_pdg(:, :, kk), corr_pdg(:, :, 1))));
% end



%% pooling data across fields

% CC bs
% initialize these cell arrays once
nat_btwses_cc_cat = cell(1, 6);
pdg_btwses_cc_cat = cell(1, 6);

% run this for each mouse after calculating its values above
weeks = [1 1 1 1 1 1];                  % starting on week 2, indicates weeks for which this mouse has data. e.g. [1 1 1 1 0 0] means mouse was recorded for only 5 weeks.
if ~isempty('weeks')
    num_sessions = length(weeks);
    ct = 1;
    for kk = 1:num_sessions
        if weeks(kk)
            nat_btwses_cc_cat{kk} = [nat_btwses_cc_cat{kk} nat_btwses_cc(ct)];
            pdg_btwses_cc_cat{kk} = [pdg_btwses_cc_cat{kk} pdg_btwses_cc(ct)];
            ct = ct + 1;
        end
    end
end

% save data once all mice have been added
save nat_btwses_cc_cat nat_btwses_cc_cat
save pdg_btwses_cc_cat pdg_btwses_cc_cat


%% plotting pooled data

% C
pdgabs = importdata('pdg_meanconn_absdiff_cat.mat');
natabs = importdata('nat_meanconn_absdiff_cat.mat');

figure
hold on
for kk = 1:length(pdgabs)
    plot(1:2, [pdgabs(kk) natabs(kk)], 'k');
    scatter(1, pdgabs(kk), 'filled', 'MarkerFaceColor', 'b');
    scatter(2, natabs(kk), 'filled', 'MarkerFaceColor', 'r');
end
[~, p, ~, stats] = ttest(pdgabs, natabs);
xlabel('Stimulus')
ylabel('Mean |D42 mean conn. - D0 mean conn.|')
title(sprintf('p = %.2e', p))
ylim([0 0.12])
xlim([0 3])
set(gcf, 'Position', [200 200 250 500])

% D
nat = importdata('nat_btwses_cc_cat.mat');
pdg = importdata('pdg_btwses_cc_cat.mat');


nat = cellfun(@(x)1-x, nat, 'UniformOutput', false);
pdg = cellfun(@(x)1-x, pdg, 'UniformOutput', false);
mean_nat = cellfun(@mean, nat);
mean_pdg = cellfun(@mean, pdg);
SEM_nat = cellfun(@std, nat)./sqrt(cellfun(@length, nat));
SEM_pdg = cellfun(@std, pdg)./sqrt(cellfun(@length, pdg));

for kk = 1:length(nat)
    [~, p(kk), ~, stats(kk)] = ttest(nat{kk}, pdg{kk});
end

figure
hold on
errorbar(1:6, mean_pdg, SEM_pdg, 'b');
errorbar(1:6, mean_nat, SEM_nat, 'r');
xlim([0 7])
ylim([0.5 1])
axis square
title(sprintf('p = %.2e ', p))






