%% Draw boundaries for each field, save data

[LayerVec, boundary234, boundary45, boundary56, xaxis_ROIdensity, yaxis_ROIdensity] = drawLayerBoundaries_V2(0.838);               %input is resolution of recording in um/pixel depending on field, 0.908 for MJG003, STDM077, STDM076 (1.2x), 0.838 for MJG013, 
save LayerVec LayerVec
save boundary234 boundary234
save boundary45 boundary45
save boundary56 boundary56
save xaxis_ROIdensity xaxis_ROIdensity
save yaxis_ROIdensity yaxis_ROIdensity

%% example field comparing layers
MJG013w1 = StabilityAnalyzer();         % initialize Analyzer
MJG013w1.importData();                  % import data ('RespData', 'RoiINFO', 'StabilityData')
MJG013w1.setQualityThreshold(3);        % set quality threshold
MJG013w1.cellSelection;                 % select inclusion criteria ('quality', 'presence', 'PDG_Responsive_thresh', 'NatMov_Responsive_thresh')
num_sessions = MJG013w1.num_sessions;
num_cells = MJG013w1.num_cells;
qual = MJG013w1.getUse_cells;
sessions = MJG013w1.getUse_sessions;
present = sum(sessions, 2) == MJG013w1.num_sessions;
good_cells = qual & present';

load('LayerVec.mat')                    % load in field's layer information
L23 = strcmp(LayerVec, 'L2/3');
L4 = strcmp(LayerVec, 'L4');
L5 = strcmp(LayerVec, 'L5');
L6 = strcmp(LayerVec, 'L6');

PDG_RDI = MJG013w1.StabilityData.PDG.RDI;
NatMov_RDI = MJG013w1.StabilityData.NatMov.RDI;

mean_PDG_RDI = mean(PDG_RDI, 1);
mean_NatMov_RDI = mean(NatMov_RDI, 1);

NatMov_RDI_L23 = mean_NatMov_RDI(good_cells & L23);
NatMov_RDI_L4 = mean_NatMov_RDI(good_cells & L4);
NatMov_RDI_L5 = mean_NatMov_RDI(good_cells & L5);

PDG_RDI_L23 = mean_PDG_RDI(good_cells & L23);
PDG_RDI_L4 = mean_PDG_RDI(good_cells & L4);
PDG_RDI_L5 = mean_PDG_RDI(good_cells & L5);

figure
boxplot([NatMov_RDI_L23 NatMov_RDI_L4 NatMov_RDI_L5], [zeros(1, length(NatMov_RDI_L23)) ones(1, length(NatMov_RDI_L4)) 2*ones(1, length(NatMov_RDI_L5))]);
[~, ~, stats] = anova1([NatMov_RDI_L23 NatMov_RDI_L4 NatMov_RDI_L5], [zeros(1, length(NatMov_RDI_L23)) ones(1, length(NatMov_RDI_L4)) 2*ones(1, length(NatMov_RDI_L5))], 'off');
c = multcompare(stats);

figure
boxplot([PDG_RDI_L23 PDG_RDI_L4 PDG_RDI_L5], [zeros(1, length(PDG_RDI_L23)) ones(1, length(PDG_RDI_L4)) 2*ones(1, length(PDG_RDI_L5))]);
[~, ~, stats] = anova1([PDG_RDI_L23 PDG_RDI_L4 PDG_RDI_L5], [zeros(1, length(PDG_RDI_L23)) ones(1, length(PDG_RDI_L4)) 2*ones(1, length(PDG_RDI_L5))], 'off');
c = multcompare(stats);

%% ROI patches pseudocolored by average RDI
load NatMov_ChronicImaging_maps.mat;     % load in field maps for example field
load('TSeries-01172019-1737-W1-NatMov-002_registered_warped_data.mat')         % load in example (unprocessed) 2P data file to reference cell masks

MJG013w1.cellSelection;               % select quality, presence, and NatMov_Responsive_thresh
num_cells = MJG013w1.num_cells;
qual = MJG013w1.getUse_cells;
sessions = MJG013w1.getUse_sessions;
present = sum(sessions, 2) == MJG013w1.num_sessions;
good_cells = qual & present';

figure
imagesc(NatMov_ChronicImaging_maps.avg_projections(:, :, 6)); colormap gray; axis square
hold on
minRDI = min(mean_NatMov_RDI(good_cells & ~L6));
mean_NatMov_RDI_scaled = mean_NatMov_RDI - minRDI;
% maxRDI = max(mean_NatMov_RDI_scaled(good_cells & ~L6));
maxRDI = 0.5;                    % or force it
mean_NatMov_RDI_scaled = mean_NatMov_RDI_scaled/maxRDI;
mean_NatMov_RDI_scaled(mean_NatMov_RDI_scaled > 1) = 1;
colors = jet;
forcemaxaxis = 'yes';
for ii = 1:length(data.cellMasks)
    if good_cells(ii) && ~strcmp(LayerVec{ii}, 'L6') && mean_NatMov_RDI(ii) >= 0
        currROI = data.cellMasks{ii};
        coloridx = ceil(mean_NatMov_RDI_scaled(ii)*64);
        if coloridx == 10
            coloridx = 1;
        end
        pdata = patch(currROI(:,1),currROI(:,2), colors(coloridx, :));
        pdata.EdgeAlpha = 0;
        pdata.FaceAlpha = .5;
    end
end
plot(boundary234(:, 1), boundary234(:, 2));
plot(boundary45(:, 1), boundary45(:, 2));
plot(boundary56(:, 1), boundary56(:, 2));
set(gca, 'xtick', [])
set(gca, 'ytick', [])

%% ROI patches pseudocolored by layer
load('NatMov_ChronicImaging_maps.mat')
load('TSeries-01172019-1737-W1-NatMov-002_registered_warped_data.mat')  % MJG013w1 

avgproj = NatMov_ChronicImaging_maps.avg_projections(:, :, 1);

MJG013w1.cellSelection;               % select quality and presence here
num_cells = MJG013w1.num_cells;
qual = MJG013w1.getUse_cells;
sessions = MJG013w1.getUse_sessions;
present = sum(sessions, 2) == MJG013w1.num_sessions;
good_cells = qual & present';

figure
imagesc(avgproj); colormap gray; axis square
hold on
for ii = 1:num_cells
    if ~strcmp(LayerVec{ii}, 'L6')
        if strcmp(LayerVec{ii}, 'L2/3')
            color = [1 0 0];
        elseif strcmp(LayerVec{ii}, 'L4')
            color = [0 1 0];
        elseif strcmp(LayerVec{ii}, 'L5')
            color = [0 0 1];
        end
        currROI = data.cellMasks{ii};
        pdata = patch(currROI(:,1),currROI(:,2), color);
        pdata.EdgeAlpha = 0;
        pdata.FaceAlpha = .5;
    end
end

set(gca, 'xtick', [])
set(gca, 'ytick', [])

%% pulling out single cell map layer examples

MJG013w1c = ChronicImaging_Mouse;
MJG013w1c.importData();
MJG013w1c.evaluateROIs;             % run 'Visual_Inpsection_A' to process cell images

num_cells = MJG013w1c.num_cells;
num_sessions = MJG013w1c.num_sessions;
neuron = 1;

box_size = size(MJG013w1c.Vis_cellInspection_maps.NatMov.cell_avgproj_mat, 1);
NatMov_maps = zeros(box_size, box_size, num_sessions, num_cells);                   % avg projections for natmov and pdg (act maps remain unaltered
PDG_maps = zeros(box_size, box_size, num_sessions, num_cells);
for ii = neuron:num_cells
    ii
    for kk = 1:num_sessions
        currmap = MJG013w1c.Vis_cellInspection_maps.NatMov.cell_avgproj_mat(:, :, kk, ii);
        currmap = imresize(currmap, [1000 1000]);
        currmap = imrotate(currmap, -90);
        currmap = imresize(currmap, [36 36]);
        NatMov_maps(:, :, kk, ii) = currmap;
        
        currmap = MJG013w1c.Vis_cellInspection_maps.PDG.cell_avgproj_mat(:, :, kk, ii);
        currmap = imresize(currmap, [1000 1000]);
        currmap = imrotate(currmap, -90);
        currmap = imresize(currmap, [36 36]);
        PDG_maps(:, :, kk, ii) = currmap;
    end
end      


neuron = 28;
map_type = 'PDG';
switch map_type
    case 'NatMov'
        maps = NatMov_maps;
    case 'PDG'
        maps = PDG_maps;
end
figure
for kk = 1:num_sessions
    subplot(1, 7, kk)
    imagesc(imadjustn(maps(:, :, kk, neuron)));
    colormap gray
    axis square
end
title(sprintf('neuron %d', neuron))



%% pooled prism fields comparing layers

% LayerVec_cat = [];            % initialize once
load('LayerVec.mat')                            % load in and add layer info for each field in same order as StabilityAnalyzer will be constructed 
LayerVec_cat = [LayerVec_cat LayerVec];
save LayerVec_cat LayerVec_cat

load('LayerVec_cat.mat')
L23 = strcmp(LayerVec_cat, 'L2/3');
L4 = strcmp(LayerVec_cat, 'L4');
L5 = strcmp(LayerVec_cat, 'L5');

allp = StabilityAnalyzer();
allp.importData();               % import first field
allp.addData();                  % add subsequent fields
allp.setQualityThreshold(3);
allp.cellSelection;             % select quality, presence, and dual responsive
num_sessions = allp.num_sessions;
num_cells = allp.num_cells;
qual = allp.getUse_cells;
sessions = allp.getUse_sessions;
sessions(isnan(sessions)) = 1;          % convert NaNs to 1s for logicals
present = sum(sessions, 2) == allp.num_sessions;
good_cells = qual;

PDG_RDI = allp.StabilityData.PDG.RDI;
NatMov_RDI = allp.StabilityData.NatMov.RDI;

mean_PDG_RDI = zeros(1, num_cells);
mean_Natmov_RDI = zeros(1, num_cells);
for ii = 1:num_cells
    mean_PDG_RDI(ii) = nanmean(PDG_RDI(logical(sessions(ii, 2:end)), ii), 1);                 % averaging across all sessions for which the cell is present
    mean_NatMov_RDI(ii) = nanmean(NatMov_RDI(logical(sessions(ii, 2:end)), ii), 1);
end

NatMov_RDI_L23 = mean_NatMov_RDI(good_cells & L23);
NatMov_RDI_L4 = mean_NatMov_RDI(good_cells & L4);
NatMov_RDI_L5 = mean_NatMov_RDI(good_cells & L5);

PDG_RDI_L23 = mean_PDG_RDI(good_cells & L23);
PDG_RDI_L4 = mean_PDG_RDI(good_cells & L4);
PDG_RDI_L5 = mean_PDG_RDI(good_cells & L5);

figure
boxplot([NatMov_RDI_L23 NatMov_RDI_L4 NatMov_RDI_L5], [zeros(1, length(NatMov_RDI_L23)) ones(1, length(NatMov_RDI_L4)) 2*ones(1, length(NatMov_RDI_L5))], 'Symbol', 'o');
ylim([-0.2 1])
axis square

figure
boxplot([PDG_RDI_L23 PDG_RDI_L4 PDG_RDI_L5], [zeros(1, length(PDG_RDI_L23)) ones(1, length(PDG_RDI_L4)) 2*ones(1, length(PDG_RDI_L5))], 'Symbol', 'o');
ylim([-0.2 1])
axis square


% 2 way anova for testing across layers and stimuli
y = [PDG_RDI_L23 PDG_RDI_L4 PDG_RDI_L5 NatMov_RDI_L23 NatMov_RDI_L4 NatMov_RDI_L5]';
clear cell
n = length(PDG_RDI_L23) + length(PDG_RDI_L4) + length(PDG_RDI_L5);
S = cell(n*2, 1);
S(1:n) = {'PDG'}; S(n+1:end) = {'NatMov'};
La = cell(1, length(PDG_RDI_L23)); La(:) = {'L23'}; 
Lb = cell(1, length(PDG_RDI_L4)); Lb(:) = {'L4'};
Lc = cell(1, length(PDG_RDI_L5)); Lc(:) = {'L5'};
L = [La Lb Lc La Lb Lc]';

[p, tbl, stats] = anovan(y, {L S});

figure
c = multcompare(stats);
