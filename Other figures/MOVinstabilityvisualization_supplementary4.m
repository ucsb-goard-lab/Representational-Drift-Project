% NatMov pixeltuning visualization

%generating the maps
colorimg = NaturalMovie_PixelTuning;

imagesc(data.colorimg);
axis square; box off; set(gca, 'yTick', []); set(gca, 'xTick', []); 




% visualize ROI IDs for selecting example cells
load('TSeries-04042019-1442-NatMov-003_registered_warped_data.mat');

load('RoiINFO.mat');        %import cell quality info for this field

figure;
imagesc(data.colorimg);
axis square       
for ii = 1:size(data.DFF, 1)
    if RoiINFO.quality(ii) >= 3 && sum(RoiINFO.presence(ii, :)) == size(RoiINFO.presence, 2) && RoiINFO.NatMov_Responsive_thresh(ii)
%         pdata = patch(data(1).cellMasks{1,ii}(:, 1), data(1).cellMasks{1,ii}(:, 2), 'c');
%         pdata.EdgeAlpha = 0;
%         pdata.FaceAlpha = .5;
        cellID = text(data(1).cellMasks{1,ii}(10, 1), data(1).cellMasks{1,ii}(10, 2), num2str(ii), 'Color', 'w', 'Fontsize', 9);
        hold on
    end
end

% example cell RDI curves
X147A = ResponseVisualizer;
X147A.importData();
X147A.setQualityThreshold(3);
X147A.cellSelection;
qual = X147A.getUse_cells;
sessions = X147A.getUse_sessions;
present = sum(sessions, 2) == X147A.num_sessions;
good_cells = qual' & present;

neuron = 692;     
figure
for ii = neuron:X147A.num_cells
    if good_cells(ii)
        NatMovcurve = X147A.StabilityData.NatMov.RDI(:, ii);
        plot(NatMovcurve, 'color', [0.8 0.5 0.5])
        ylim([-0.2 1])
        xlim([0 X147A.num_sessions])
        refline(0, X147A.StabilityData.NatMov.RDI_control(ii))
        axis square; box off
        title(sprintf('cell %d', ii))
        pause
    end  
end


% example cell average traces
X147A.plotPDGvNatMov_avgs(neuron, [0 0.5]);
set(gcf, 'Position', [400 50 800 1200])

%extract trial-averages for example neuron
NatMov_mean = squeeze(mean(X147A.RespData.NatMov.RespMat_onTime(:, :, neuron, :), 1));