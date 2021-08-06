classdef ResponseVisualizer < Analyzer
    
    properties
        normalized_responses
    end
    
    methods
        
        function obj = ResponseVisualizer()
        end
        
        function obj = plotPDG_traces(obj, startcell)  
            %input cells to iterate through (useCells) and cell to start with (startcell)
            %Defaults are to use cells responsive to gratings and start on cell 1
            %see plotSingleStimulus_traces method
            if nargin < 2
                startcell = 1;
            end

            obj.plotSingleStimulus_traces('PDG', startcell)
        end

        function obj = plotNatMov_traces(obj, startcell)
            %input cells to iterate through (useCells) and cell to start with (startcell)
            %Defaults are to use cells responsive to movies and start on cell 1
            %see plotSingleStimulus_traces method
            if nargin < 3
                startcell = 1;
            end

            obj.plotSingleStimulus_traces('NatMov', startcell)
        end

        function obj = plotPDGvNatMov_traces2(obj, startcell, norm_flag)
            %Displays PDG trial to trial activity alongside NatMov activity, and comparison of NatMov average traces (far right)
            %Input cells to iterate through, and cell to start on
            %Defaults are cells responsive to both gratings and movies (lenient) and cell 1
            if nargin < 2
                startcell = 1;
            end
            
            if isempty(obj.quality_threshold)
                fprintf('ERROR: Must set quality threhold first.\n');
                return
            end

            gain = 80;

            PDG_RespMat = obj.RespData.PDG.RespMat_Full;
            % PDG_CCs = obj.StabilityData.PDG.CCs;
            NatMov_RespMat = obj.RespData.NatMov.RespMat_onTime;
            % NatMov_CCs = obj.StabilityData.NatMov.CCs;
            combined(1).RespMat = PDG_RespMat;
            combined(2).RespMat = NatMov_RespMat;
            
            normalized = obj.normalizeALL(combined, norm_flag);
            PDG_RespMat = normalized(1).RespMat;
            NatMov_RespMat = normalized(2).RespMat;

            disp_flag = questdlg('Plot:', 'Dialog', 'All', 'Included only', 'All');

            disp('Press any key to move to next cell, click to save as png file');

            PDG_RespMat = obj.catSessions(PDG_RespMat);
            NatMov_RespMat = obj.catSessions(NatMov_RespMat);

            flags = obj.getUse_cells;
            sessions = obj.getUse_sessions;

            for ii = startcell:obj.num_cells
                if strcmp(disp_flag, 'All') || flags(ii)
                    f = figure;
                    subplot(1, 3, 1)
                    image(PDG_RespMat(:, :, ii)*gain);
                    title(sprintf('PDG neuron %d', ii));
                    xticks([])
                    yticks([])
                    xlabel('Frames')
                    ylabel('Trials')
                    box off

                    subplot(1, 3, 2)
                    image(NatMov_RespMat(:, :, ii)*gain);
                    titlestring = sprintf('%d, ', flags(ii));
                    for kk = 1:obj.num_sessions
                        titlestring = sprintf('%s %d', titlestring, sessions(ii, kk));
                    end
                    title(sprintf('NatMov neuron %d, %s', ii, titlestring));
                    xticks([])
                    yticks([])
                    xlabel('Frames')
                    ylabel('Trials')
                    box off

                    subplot(1, 3, 3)
                    if obj.num_sessions > 2
                        plot(obj.StabilityData.PDG.RDI(:, ii), 'Color', 'b')
                        hold on
                        plot(obj.StabilityData.NatMov.RDI(:, ii), 'Color', 'r');
                        xlim([0 obj.num_sessions+1])
                        ylim([0 1])
                    else
                        scatter(1, obj.StabilityData.PDG.RDI(:, ii), 60, 'filled', 'MarkerFaceColor', 'b');
                        hold on
                        scatter(1, obj.StabilityData.NatMov.RDI(:, ii), 60, 'filled', 'MarkerFaceColor', 'r');
                        if obj.StabilityData.NatMov.RDI(:, ii) < 0 || obj.StabilityData.PDG.RDI(:, ii) < 0
                            ylim([-1 1])
                        else
                            ylim([0 1])
                        end
                    end


                    set(f,'color',[1 1 1])
                    set(f,'Position',[400 20 1800 800])
                    w = waitforbuttonpress;
                    if ~w
                        saveas(gcf, ['n' num2str(ii) '.png']);
                        uistack(f, 'top')
                    end 
                    close
                end
            end
        end

        function obj = plotTempoStructure_traces(obj, startcell, norm_flag)
            %Displays PDG trial to trial activity alongside NatMov activity, and comparison of NatMov average traces (far right)
            %Input cells to iterate through, and cell to start on
            %Defaults are cells responsive to both gratings and movies (lenient) and cell 1
            if nargin < 2
                startcell = 1;
            end
            
            if isempty(obj.quality_threshold)
                fprintf('ERROR: Must set quality threhold first.\n');
                return
            end
            
            
            gain = 80;

            PDG_RespMat = obj.RespData.PDG.RespMat_Full;
            % PDG_CCs = obj.StabilityData.PDG.CCs;
            NatMov_RespMat = obj.RespData.NatMov.RespMat_onTime;
            % NatMov_CCs = obj.StabilityData.NatMov.CCs;
            PDGcont_RespMat = obj.RespData.PDGcont.RespMat_onTime;
            NatMovdisc_RespMat = obj.RespData.NatMovdisc.RespMat_Full;
            combined(1).RespMat = PDG_RespMat;
            combined(2).RespMat = NatMov_RespMat;
            combined(3).RespMat = PDGcont_RespMat;
            combined(4).RespMat = NatMovdisc_RespMat;

            normalized = obj.normalizeALL(combined, norm_flag);
            PDG_RespMat = normalized(1).RespMat;
            NatMov_RespMat = normalized(2).RespMat;
            PDGcont_RespMat = normalized(3).RespMat;
            NatMovdisc_RespMat = normalized(4).RespMat;

            disp_flag = questdlg('Plot:', 'Dialog', 'All', 'Included only', 'All');

            disp('Press any key to move to next cell, click to save as png file');

            PDG_RespMat = obj.catSessions(PDG_RespMat);
            NatMov_RespMat = obj.catSessions(NatMov_RespMat);
            PDGcont_RespMat = obj.catSessions(PDGcont_RespMat);
            NatMovdisc_RespMat = obj.catSessions(NatMovdisc_RespMat);


            flags = obj.getUse_cells;
            sessions = obj.getUse_sessions;

            for ii = startcell:obj.num_cells
                if strcmp(disp_flag, 'All') || flags(ii)
                    f = figure;
                    subplot(1, 4, 1)
                    image(PDG_RespMat(:, :, ii)*gain);
                    title(sprintf('PDG neuron %d', ii));
                    xticks([])
                    yticks([])
                    xlabel('Frames')
                    ylabel('Trials')
                    box off

                    subplot(1, 4, 2)
                    image(NatMov_RespMat(:, :, ii)*gain);
                    titlestring = sprintf('%d, ', flags(ii));
                    for kk = 1:obj.num_sessions
                        titlestring = sprintf('%s %d', titlestring, sessions(ii, kk));
                    end
                    title(sprintf('NatMov neuron %d, %s', ii, titlestring));
                    xticks([])
                    yticks([])
                    xlabel('Frames')
                    ylabel('Trials')
                    box off

                    subplot(1, 4, 3)
                    image(PDGcont_RespMat(:, :, ii)*gain);
                    titlestring = sprintf('%d, ', flags(ii));
                    for kk = 1:obj.num_sessions
                        titlestring = sprintf('%s %d', titlestring, sessions(ii, kk));
                    end
                    title(sprintf('PDGcont neuron %d, %s', ii, titlestring));
                    xticks([])
                    yticks([])
                    xlabel('Frames')
                    ylabel('Trials')
                    box off

                    subplot(1, 4, 4)
                    image(NatMovdisc_RespMat(:, :, ii)*gain);
                    titlestring = sprintf('%d, ', flags(ii));
                    for kk = 1:obj.num_sessions
                        titlestring = sprintf('%s %d', titlestring, sessions(ii, kk));
                    end
                    title(sprintf('NatMovdisc neuron %d, %s', ii, titlestring));
                    xticks([])
                    yticks([])
                    xlabel('Frames')
                    ylabel('Trials')
                    box off


                    set(f,'color',[1 1 1])
                    set(f,'Position',[400 20 1800 800])
                    w = waitforbuttonpress;
                    if ~w
                        saveas(gcf, ['n' num2str(ii) '.png']);
                        uistack(f, 'top')
                    end 
                    close
                end
            end
        end

        function obj = plotPDGvNatMov_avgs(obj, startcell, ylims)
            if nargin < 2
                startcell = 1;
            end

            if nargin < 3
                ylims = [];
            end

            PDG_RespMat = obj.RespData.PDG.RespMat_Full;
            PDG_CCs = obj.StabilityData.PDG.CCs;
            NatMov_RespMat = obj.RespData.NatMov.RespMat_onTime;
            NatMov_CCs = obj.StabilityData.NatMov.CCs;

            [PDG_RespMat, NatMov_RespMat] = obj.normalizeALL(PDG_RespMat, NatMov_RespMat, 'minimum');

            PDG_onFrames = size(PDG_RespMat, 2);
            NatMov_onFrames = size(NatMov_RespMat, 2);

            flags = obj.getUse_cells;
            sessions = obj.getUse_sessions;

            figure
            for ii = startcell:obj.num_cells
                PDG_mean = squeeze(mean(PDG_RespMat(:, :, ii, :), 1));
                NatMov_mean = squeeze(mean(NatMov_RespMat(:, :, ii, :), 1));

                for kk = 1:obj.num_sessions
                    
                    subplot(obj.num_sessions, 2, kk*2-1)
                    plot(PDG_mean(:, kk));
                    if ~isempty(ylims)
                        ylim([ylims(1) ylims(2)])
                    else
                        ylim([min(min(PDG_mean)) max(max(PDG_mean))]);
                    end
                    title(sprintf('Neuron %d, %d', ii, flags(ii)));

                    subplot(obj.num_sessions, 2, kk*2)
                    plot(NatMov_mean(:, kk));
                    if ~isempty(ylims)
                        ylim([ylims(1) ylims(2)])
                    else
                        ylim([min(min(PDG_mean)) max(max(PDG_mean))]);
                    end                    
                    title(sprintf('%d', sessions(ii, kk)));

                end
                pause
            end
        end
        
        function obj = plotSingleStimulus_traces(obj, stimulus_type, startcell, norm_flag)               %this method is shared between PDG and NatMov trace plotting methods
            %displays response matricies alongside average traces of each session, reference session overlayed with subsequent session
            switch stimulus_type
                case 'PDG'
                    RespMat = obj.RespData.PDG.RespMat_Full;
                    onFrames = size(RespMat, 2);
                    CCs = obj.StabilityData.PDG.CCs;
                    try cd('PDG_traces')
                    catch
                        mkdir('PDG_traces')
                        cd('PDG_traces')
                    end
                case 'NatMov'
                    RespMat = obj.RespData.NatMov.RespMat_onTime;
                    onFrames = size(RespMat, 2);
                    CCs = obj.StabilityData.NatMov.CCs;
                    try cd('NatMov_traces')
                    catch
                        mkdir('NatMov_traces')
                        cd('NatMov_traces')
                    end
            end

            if nargin < 4
                norm_flag = 'no';
            end

            RespMat = obj.catSessions(RespMat);

            if ~strcmp(norm_flag, 'no')
                [~, RespMat] = obj.normalizeALL(RespMat, RespMat, norm_flag);
            end

            screening = questdlg('Do you wish to manually curate the cells?', 'Query', 'Curate', 'Autosave', 'Curate');

            if strcmp(screening, 'Curate')
                disp('Press any key to move to next cell, click to save as png file');
            else
                disp('Autosaving all figures... see you in a bit. Go take a walk or something.')
            end
            
            obj.cellSelection;
            cell_vec = obj.getUse_cells;
            sessions = obj.getUse_sessions;
            
            figure;
            for ii = startcell:obj.num_cells
                if cell_vec(ii)
                    imagesc(RespMat(:, :, ii));
                    titlestring = sprintf('%d, ', cell_vec(ii));
                    for kk = 1:obj.num_sessions
                        titlestring = sprintf('%s %d', titlestring, sessions(ii, kk));
                    end
                    title(sprintf('NatMov neuron %d, %s', ii, titlestring));                    
                    ylabel('Trials')
                    set(gcf,'color',[1 1 1])
                    set(gcf,'Position',[20 20 1000 1400])
                    set(gca, 'xtick', []);
                    set(gca, 'ytick', []);
                    box off
                    if strcmp(screening, 'Curate')
                        w = waitforbuttonpress;
                        if ~w
                            saveas(gcf, sprintf('n%d.png', ii))
                        end
                    else
                        saveas(gcf, sprintf('n%d.png', ii))
                        pause(1)
                    end   
                end
            end
            if strcmp(screening, 'Curate')
                fprintf('\n ...Welcome back. %d figures saved.\n', sum(cell_vec))
                close
            end
        end
            
        function outRespMats = normalizeALL(obj, inRespMats, method)
            if ~isempty(obj.normalized_responses)
                overwrite = questdlg('Overwrite?', 'Normalized data detected, overwrite?', 'Yes', 'No', 'Yes');
            else
                overwrite = 'Yes';
            end

            if strcmp(overwrite, 'Yes');
                outRespMats = inRespMats;
                num_Mats = length(inRespMats);
                for ii = 1:obj.num_cells
                    switch method
                    case 'minimum'
                        mins = zeros(1, num_Mats);
                        for rr = 1:num_Mats
                            mins(rr) = min(min(min(squeeze(inRespMats(rr).RespMat(:, :, ii, :)))));
                        end
                        MIN = min(mins);
                    case 'first mode'
                        modes = zeros(1, num_Mats);
                        for rr = 1:num_Mats
                            modes(rr) = mode(squeeze(inRespMats(rr).RespMat(:, :, ii, :)));
                        end
                        MIN = min(modes);
                    end

                    for rr = 1:num_Mats
                        inRespMats(rr).RespMat(:, :, ii, :) = inRespMats(rr).RespMat(:, :, ii, :) - MIN;
                    end
                    maxes = zeros(1, num_Mats);
                    for rr = 1:num_Mats
                        maxes(rr) = max(max(max(squeeze(inRespMats(rr).RespMat(:, :, ii, :)))));
                    end
                    MAX = max(maxes);

                    for rr = 1:num_Mats
                        outRespMats(rr).RespMat(:, :, ii, :) = inRespMats(rr).RespMat(:, :, ii, :)/MAX;
                    end
                    ii
                end
                obj.normalized_responses = outRespMats;
            else
                outRespMats = obj.normalized_responses;
            end
        end

        function peakResponseColoring(obj, RespMat, session)
            qual = obj.getUse_cells;
            sessions = obj.getUse_sessions;
            good_cells = qual & sessions';

            fprintf('Select datafile\n');
            filename = uigetfile('.mat');
            datafile = importdata(filename);

            fprintf('Select maps file\n');
            filename = uigetfile('.mat');
            maps = importdata(filename);

            cellMasks = datafile.cellMasks;
            responsive = obj.RoiINFO.NatMov_Responsive_shuffle & obj.RoiINFO.NatMov_Responsive_shuffle;

            mean_RespMat = squeeze(mean(RespMat, 1));
            peak_idx = zeros(1, obj.num_cells);
            for ii = 1:obj.num_cells
                [~, peak_frame] = max(mean_RespMat(:, ii));
                peak_idx(ii) = peak_frame/size(mean_RespMat, 1);     %normalized to the number of frames
            end

            imagesc(maps.avg_projections(:, :, session)); axis square; colormap gray;
            hold on
            colors = jet;
            for ii = 1:obj.num_cells
                if good_cells(session, ii) && responsive(ii)
                    currROI = cellMasks{ii};
                    coloridx = ceil(peak_idx(ii)*64);
                    if coloridx == 0
                        coloridx = 1;
                    end
                    pdata = patch(currROI(:,1),currROI(:,2), colors(coloridx, :));
                    pdata.EdgeAlpha = 0;
                    pdata.FaceAlpha = .5;
                end
            end
            set(gcf,'color',[1 1 1])
            set(gcf,'Position',[300 100 1000 1000])
            box off
            set(gca, 'xtick', []);
            set(gca, 'ytick', []);
        end

        function out = catSessions(obj, in)
            out = [];
            for kk = 1:obj.num_sessions
                out = cat(1, out, in(:, :, :, kk));
            end
        end
    end
end