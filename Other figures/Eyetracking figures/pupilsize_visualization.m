


%% pupil size and gain effect visualization

% MJG033 for frame examples
eyemovie = EyeTracker();
eyemovie.readEyeTrackingVideo();
eyemovie.cleanVideo('interpolate');
movie = eyemovie.movie;

start = 3000;
figure; colormap gray; box off
for tt = start:15000
    imagesc(movie(:, :, tt));
    title(sprintf('frame %d', tt))
%     pause
    pause(0.01)
end

example_large = movie(60:150, 90:190, 200);
figure; colormap gray; box off;
imagesc(example_large);
    
example_small = movie(60:150, 90:190, 3180);
figure; colormap gray; box off;
imagesc(example_small);

load('RespData.mat');
load('pupilinfo_MOV.mat');

% sort trials by pupil area
area = processed_pupil.area(1).sorted;
mean_area = mean(area(:, 51:350), 2);
[~, idx] = sort(mean_area, 'descend');

MOV = RespData.NatMov.RespMat_onTime(:, :, :, 1);
% import a mouse's data
vis = ResponseVisualizer();
vis.importData();
vis.setQualityThreshold(3);
vis.cellSelection;
qual = vis.getUse_cells;

bigpupil_avg = squeeze(mean(MOV(idx(1:10), :, :), 1));
bigpupil_sem = squeeze(std(MOV(idx(1:10), :, :), [], 1))/sqrt(10);
smallpupil_avg = squeeze(mean(MOV(idx(21:end), :, :), 1));
smallpupil_sem = squeeze(std(MOV(idx(21:end), :, :), [], 1))/sqrt(10);

for ii = 1:length(qual)
    if qual(ii)
        figure
        subplot(2, 1, 1)
        hold on
        plot(bigpupil_avg(:, ii))
        plot(bigpupil_avg(:, ii) - bigpupil_sem(:, ii))
        plot(bigpupil_avg(:, ii) + bigpupil_sem(:, ii))
        hold off
        ylim([min([min(bigpupil_avg(:, ii) - bigpupil_sem(:, ii)) min(smallpupil_avg(:, ii) - smallpupil_sem(:, ii))]) max([max(bigpupil_avg(:, ii) + bigpupil_sem(:, ii)) max(smallpupil_avg(:, ii) + smallpupil_sem(:, ii))])]) 
        title(sprintf('neuron %d', ii))
        subplot(2, 1, 2)
        hold on
        plot(smallpupil_avg(:, ii))
        plot(smallpupil_avg(:, ii) - smallpupil_sem(:, ii))
        plot(smallpupil_avg(:, ii) + smallpupil_sem(:, ii))
        hold off
        ylim([min([min(bigpupil_avg(:, ii) - bigpupil_sem(:, ii)) min(smallpupil_avg(:, ii) - smallpupil_sem(:, ii))]) max([max(bigpupil_avg(:, ii) + bigpupil_sem(:, ii)) max(smallpupil_avg(:, ii) + smallpupil_sem(:, ii))])]) 
        pause
        close
    end
end
       