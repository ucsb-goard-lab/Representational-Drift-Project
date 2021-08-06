% gratings stability supplementary

%% Pixeltuning visualization
raw_merge = OrientationTuning_PixelTuning;
save raw_merge_SKKS058_D42 raw_merge

new = imgaussfilt(raw_merge);
image(new)



%% OSI distributions over all sessions 

an = StabilityAnalyzer();
an.importData();
an.addData();
an.addData([1 1 0 1 1 1]);            % vector indicating weeks for which data is present for the given mouse. [1 1 1 0 1 1 1] indicates data corresponds to weeks 1, 2, 3, 5, 6, 7.
an.setQualityThreshold(3);
an.cellSelection;               % select 'quality' and 'presence'
num_sessions = an.num_sessions;
num_cells = an.num_cells;
qual = an.getUse_cells;          
sessions = an.getUse_sessions;
present = sum(sessions, 2) == an.num_sessions;
isTuned = an.OrientationData.isTuned;
good_cells = qual & present' & isTuned(1, :);           % quality, present on all sessions, and orientation tuned for first session

osiMat = an.OrientationData.osiMat;
osiMat = osiMat(:, good_cells)';
% violinplot(osiMat, [], 'ShowViolin', false, 'ShowMedian', false, 'ShowWhisker', false, 'ShowBox', false', 'ShowNotches', false, 'ViolinAlpha', 0.4);
% ylim([0 1])
% hold on

mean_osi = nanmean(osiMat, 1);
SEM_osi = nanstd(osiMat, [], 1)/sqrt(sum(good_cells));

scatter(1:num_sessions, mean_osi, 'filled');
errorbar(1:num_sessions, mean_osi, SEM_osi, 'LineStyle', 'none');


%% OSI weeks apart
% perform and save data for each mouse
an = StabilityAnalyzer();
an.importData();
an.setWeeks([1 1 1 1 1 1 1]);       % set for each mouse, see above
an.setQualityThreshold(3);
an.cellSelection;
weeks_apart = an.weeksApartOSI();
save weeks_apart_osi weeks_apart

% for pooling mice, load in each mouse's data and concatenate
%initialize once
weeks_apart_osi = cell(1, 6);

load('weeks_apart_osi.mat');
for kk = 1:length(weeks_apart)
    weeks_apart_osi{kk} = [weeks_apart_osi{kk} weeks_apart{kk}];
end

mean_bysession = cellfun(@mean, weeks_apart_osi);
sem_bysession = cellfun(@std, weeks_apart_osi)./sqrt(cellfun(@length, weeks_apart_osi));

% one mouse
figure
hold on
for kk = 1:length(weeks_apart)
    scatter(kk*ones(1, length(weeks_apart{kk})), weeks_apart{kk}, 'filled', 'MarkerFaceColor', 'b');
end
xlim([0 an.num_sessions])
ylim([0 0.1])

% pooled data
figure
hold on
for kk = 1:length(weeks_apart_osi)
    scatter(kk*ones(1, length(weeks_apart_osi{kk})), weeks_apart_osi{kk}, 'filled', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.3);
end
scatter(1:length(weeks_apart_osi), mean_bysession, 'filled', 'MarkerFaceColor', 'r');
errorbar(1:length(weeks_apart_osi), mean_bysession, sem_bysession, 'LineStyle', 'none');
xlim([0 an.num_sessions])
ylim([0 0.18])


cat = [];
for kk = 1:length(weeks_apart_osi)
    cat = [cat weeks_apart_osi{kk}];
end
for kk = 1:length(weeks_apart_osi)
    [~, p(kk), ~, stats(kk)] = ttest(weeks_apart_osi{kk}, mean(cat));
end

% finding each neuron's average change in OSI across sessions
an = StabilityAnalyzer();
an.importData();
an.addData();
an.setQualityThreshold(3);
an.cellSelection;

qual = an.getUse_cells;
sessions = an.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == an.num_sessions;
isTuned = sum(an.OrientationData.isTuned, 1) == an.num_sessions;
good_cells = qual & present' & isTuned;


osiMat = an.OrientationData.osiMat;
sesavgOSIdelta = nanmean(abs(osiMat(2:end, :) - osiMat(1, :))./osiMat(1, :), 1);
meanOSIdelta = mean(sesavgOSIdelta(good_cells));                    % mean session-averaged delta OSI across all neurons tuned on all sessions
semOSIdelta = std(sesavgOSIdelta(good_cells))/sqrt(sum(good_cells));


%% oripref comparison
antest = StabilityAnalyzer();
antest.importData();
antest.addData();
antest.addData([1 1 0 1 1 1]);      % see above. [1 1 1 0 1 1 1] indicates data corresponds to weeks 1, 2, 3, 5, 6, 7. [1 1 0 1 1 1] indicates data corresponds to weeks 1, 2, 4, 5, 6. etc.
antest.setQualityThreshold(3);
antest.cellSelection;                   % select 'quality' and 'presence'
% num_sessions = antest.num_sessions;
% num_cells = antest.num_cells;
qual = antest.getUse_cells;          
sessions = antest.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == antest.num_sessions;
isTuned = antest.OrientationData.isTuned;
good_cells = qual & present' & isTuned(1, :);    

oriResp = antest.OrientationData.oriResp;
oriResp = oriResp(:, good_cells, :);

ori1 = (squeeze(oriResp(1, :, 1:6)) + squeeze(oriResp(1, :, 7:12)))/2;     
ori2 = (squeeze(oriResp(end, :, 1:6)) + squeeze(oriResp(end, :, 7:12)))/2;

%center gaussian
[~,maxIdx] = max(ori1, [], 2);
ori1_offset = 4 - maxIdx;
[~,maxIdx] = max(ori2, [], 2);
ori2_offset = 4 - maxIdx;
ori1_centered = zeros(sum(good_cells), 6);
ori2_centered = zeros(sum(good_cells), 6);
for ii = 1:sum(good_cells)
    ori1_centered(ii, :) = circshift(ori1(ii, :), ori1_offset(ii, :));
    ori2_centered(ii, :) = circshift(ori2(ii, :), ori2_offset(ii, :));
end
ori1_centered = [ori1_centered ori1_centered(:, 1)];
ori2_centered = [ori2_centered ori2_centered(:, 1)];
 
%Upsample
upsample = 30;
ori1_upsamp = zeros(sum(good_cells), 181);
ori2_upsamp = zeros(sum(good_cells), 181);
for ii = 1:sum(good_cells)
    ori1_upsamp(ii,:) = interp1(1:7, ori1_centered(ii,:), linspace(1, 7, 181));
    ori2_upsamp(ii,:) = interp1(1:7, ori2_centered(ii,:), linspace(1, 7, 181));
end

%Fit
ori1_pref = zeros(1, sum(good_cells));
ori2_pref = zeros(1, sum(good_cells));
fail_flags = zeros(1, sum(good_cells));
for ii = 1:sum(good_cells)
    ii
    try
    %    ori1_fit = fit((0:180)', ori1_upsamp(ii,:)','gauss1');
        ori1_fit = fitgauss((0:180)', ori1_upsamp(ii,:)');
    %    ori1_pref(ii) = ori1_fit.b1 - upsample*ori1_offset(ii);
        ori1_pref(ii) = ori1_fit.Coefficients.Estimate(2) - upsample*ori1_offset(ii);
        if ori1_pref(ii) < 0
            ori1_pref(ii) = ori1_pref(ii) + 180;
        end
 
    %    ori2_fit = fit((0:180)', ori2_upsamp(ii,:)','gauss1');
        ori2_fit = fitgauss((0:180)', ori2_upsamp(ii,:)');
    %    ori2_pref(ii) = ori2_fit.b1 - upsample*ori2_offset(ii);
        ori2_pref(ii) = ori2_fit.Coefficients.Estimate(2) - upsample*ori2_offset(ii);
        if ori2_pref(ii) < 0
            ori2_pref(ii) = ori2_pref(ii) + 180;
        end
    catch
        fail_flags(ii) = 1;
    end
end
greater = ori1_pref > 180 | ori2_pref > 180;
less = ori1_pref < -180 | ori2_pref < -180;
ori1_pref(greater | less) = NaN;
ori2_pref(greater | less) = NaN;

% osi_colormap = zeros(10, 3);
% osi_colormap(:, 3) = (1:10)/10;
osiMat = antest.OrientationData.osiMat;
osiVec = osiMat(1, good_cells)';
osiVec(greater | less) = NaN;
% colormap_idx = floor(rescale(osiMat(~(greater | less), 1)*100, 1, 10));
oriscatter.ori1_pref = ori1_pref;
oriscatter.ori2_pref = ori2_pref;
oriscatter.osiVec = osiVec;
save oriscatter oriscatter

figure
scatter(ori1_pref, ori2_pref, 36, osiVec, 'filled');
axis square
colormap copper
refline(1, -30)
refline(1, 150)
refline(1, -150)
refline(1, 30)
xlim([0 180])
ylim([0 180])

%combining mice (oriscatter)
cat_oriscatter = struct;
cat_oriscatter.ori1_pref = [];
cat_oriscatter.ori2_pref = [];
cat_oriscatter.osiVec = [];
load('oriscatter.mat');
cat_oriscatter.ori1_pref = cat(2, cat_oriscatter.ori1_pref, oriscatter.ori1_pref);
cat_oriscatter.ori2_pref = cat(2, cat_oriscatter.ori2_pref, oriscatter.ori2_pref);
cat_oriscatter.osiVec = cat(1, cat_oriscatter.osiVec, oriscatter.osiVec);

save cat_oriscatter cat_oriscatter

figure
scatter(cat_oriscatter.ori1_pref, cat_oriscatter.ori2_pref, 60, cat_oriscatter.osiVec, 'filled');
axis square
colormap copper
refline(1, -30)
refline(1, 150)
refline(1, -150)
refline(1, 30)
xlim([0 180])
ylim([0 180])
colorbar
caxis([0.4 1])




% individual tuning curves
SKKS058 = StabilityAnalyzer();
SKKS058.importData();
SKKS058.setQualityThreshold(3);
SKKS058.cellSelection;

oriResp = SKKS058.OrientationData.oriResp;
num_sessions = SKKS058.num_sessions;
num_cells = SKKS058.num_cells;
qual = SKKS058.getUse_cells;          
sessions = SKKS058.getUse_sessions;
present = sum(sessions, 2) == SKKS058.num_sessions;
isTuned = SKKS058.OrientationData.isTuned;
good_cells = qual & present' & isTuned(1, :);  

for ii = 1:num_cells  
    if good_cells(ii)
        figure
        curr_oriresp = squeeze(oriResp(:, ii, :));
        theta = [0 2*pi];
        pho = [20 20]; % set min
        h1 = polarplot(theta,pho); % plot white point
        set(h1,'color',[1 1 1])
        theta = [0 2*pi];
        hold on
        num_ori = size(curr_oriresp, 2);
        upsample = 12;
        sample1 = 1:upsample:num_ori*upsample;
        sample2 = 1:1:num_ori*upsample;
        theta = [0:360/size(curr_oriresp, 2)/upsample:360]*pi/180;
        pho = zeros(length(theta), 7);
        for kk = 1:num_sessions
            curr_pho = cat(2, interp1(cat(2, sample1, upsample*num_ori), cat(2, squeeze(curr_oriresp(kk,:)), curr_oriresp(kk, 1)), sample2, 'linear', 'extrap'), curr_oriresp(kk,1))';
            curr_pho(curr_pho<0) = 0; 
            h1 = polarplot(theta,curr_pho);
            pho(:, kk) = curr_pho;
        end
        % set(h1,'linewidth',2,'color',[0 0.5 0.8])
%         rlim([0 max(max(curr_oriresp))])
        title(sprintf('neuron %d', ii))

        pause
        close
    end
end









