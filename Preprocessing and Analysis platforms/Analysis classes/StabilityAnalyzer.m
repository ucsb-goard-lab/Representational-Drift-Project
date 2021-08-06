classdef StabilityAnalyzer < Analyzer

    properties
        visualize_quality
    end

    properties (Access = protected)
        RDIdata
        n
        spikeResp                      % Deconvolved response matrix (same dimensions as normal response matrix), for use in waveform extraction/event detector
    end
    
    properties (Access = private)
        called_RDIs = {};    %last RDIs plotted
        MultiMov_called     %trackes how many times MultiMov was called to plot RDI (for updating colors)
    end
    
    methods
        
        function obj = StabilityAnalyzer()
            obj.visualize_quality = false;
        end
        
        function out = RDIplotter(obj, stim_type, plot_type, weeks)
            %plot_type describes data to be plotted, see cases below
            if isempty(obj.quality_threshold)
                fprintf('ERROR: Must set quality threshold first.\n');
                return
            end

            if nargin > 3
                obj.setWeeks(weeks);
                present_weeks = find(weeks);
                timepoints = present_weeks(2:end)-1;
            else
                timepoints = 1:obj.num_sessions-1;
            end
            
            obj.cellSelection;      %select cells to use from criteria list
            
            [curr_RDIdata, obj.n.(stim_type)] = obj.extractRDI(stim_type);

            obj.setRDIdata(structPacker(obj.getRDIdata, curr_RDIdata, stim_type));
            
            switch plot_type
                case 'Average'
                    if ~strcmp(stim_type, 'MultiMov')
                        out = obj.plotRDIavg(stim_type, timepoints);
                    else
                        obj.MultiMov_called = 1;            %initialize tracker
                        for yy = 1:length(obj.RDIdata.MultiMov.RDI_avg)
                            out = obj.plotRDIavg(stim_type, timepoints);               %if MultiMov stimulus, run this for all elements of the RDI_avg cell, have the submethod smartly update colors according to how many times its been called for MultiMov stimulus
                            hold on
                            obj.MultiMov_called = obj.MultiMov_called + 1;
                        end
                    end
                case 'Violin'
                    % obj.plotRDIviolin(stim_type);
            end
            for jj = 1:length(obj.n.(stim_type))
                fprintf('Cells used in session %d: %d\n', jj, obj.n.(stim_type)(jj));
            end
        end
       
        function out = plotRDIavg(obj, stim_type, timepoints)
            color = obj.getColorMap(stim_type);

            if obj.visualize_quality
                switch obj.getQualityThreshold;
                case 5
                    style = '-';
                case 4 
                    style = '--';
                case 3 
                    style = ':';
                case 2
                    style = '-.';
                end
            else
                style = '-';
            end
            if strcmp(stim_type, 'MultiMov')
                currplot = errorbar(timepoints, obj.RDIdata.(stim_type).RDI_avg{obj.MultiMov_called}, obj.RDIdata.(stim_type).RDI_SEM{obj.MultiMov_called}, 'LineWidth', 2, 'color', color, 'LineStyle', style);
            else
                currplot = errorbar(timepoints, obj.RDIdata.(stim_type).RDI_avg, obj.RDIdata.(stim_type).RDI_SEM, 'LineWidth', 2, 'color', color, 'LineStyle', style);
            end
            xlabel('Weeks after D0')
            ylabel('Instability Index (CCws - CCbs)/(CCws + CCbs)')
            legend(currplot, stim_type, 'Location', 'East')
            ylim([-0.2 1]);
            xlim([timepoints(1)-0.5 timepoints(end)+0.5]);
            title(sprintf('n = %d', obj.n.(stim_type)(1)))          %delete this when presence is used for criteria
            yticks([0 0.5 1]);
            axis square

            out = obj.RDIdata.(stim_type).RDI_avg;
        end
        
        function plotRDIviolin(obj, stim_type)
            color = obj.getColorMap(stim_type);
            alpha = 0.8;
            violinplot(obj.RDIdata.(stim_type).RDI_included_scatter', [], 'ShowViolin', false, 'ShowBox', false, 'ShowWhisker', false, 'ShowMedian', false, 'ViolinColor', color, 'ViolinAlpha', alpha);
            xlabel('Weeks after D0')
            ylabel('Instability Index (CCws - CCbs)/(CCws + CCbs)')
            ylim([-1 1]);
            xlim([0.5 obj.num_sessions-0.5]);
            title(sprintf('n = %d', obj.n.(stim_type)(1)))          %delete this when presence is used for criteria
            yticks([-1 -0.5 0 0.5 1]);
            axis square
        end

        function out = weeksApartOSI(obj)
            qual = obj.getUse_cells;
            sessions = obj.getUse_sessions;
            isTuned = obj.OrientationData.isTuned(1, :);

            max_distance = obj.num_sessions-1;
            out = cell(1, max_distance);
            weeks = obj.getWeeks;
            osiMat = obj.OrientationData.osiMat;            % [cells x sessions]
            for kk = 1:obj.num_sessions
                for qq = 1:obj.num_sessions
                    if (qq > kk) && weeks(qq) && weeks(kk)
                        curr_dist = qq-kk;
                        osi_diff = abs(osiMat(kk, :) - osiMat(qq, :));
                        out{curr_dist} = [out{curr_dist} nanmean(osi_diff(qual & isTuned & sessions(:, kk)' & sessions(:, qq)'))];
                    end
                end
            end
        end
        
        function [plotdata, stats] = responsivity(obj, test, stim, timespan, binsize)

            logscale = 0;       % flag for plotting in semi log scale

            qual = obj.getUse_cells;
            sessions = obj.getUse_sessions;
            sessions(isnan(sessions)) = 1;
            present = sum(sessions, 2) == obj.num_sessions;
            good_cells = qual' & present;

            switch timespan
            case 'first'
                switch stim
                case 'NatMov'
                    Resp(1, :, :, :) = obj.RespData.(stim).RespMat_onTime(:, :, :, 1);        %first session
                case 'PDG'
                    Resp(1, :, :, :) = obj.RespData.(stim).RespMat_Full(:, :, :, 1);
                end

                if strcmp(test, 'reliability')
                    reliability = obj.StabilityData.(stim).CCs(1, :);     %first session
                end

                Resp = permute(Resp, [3 4 1 2]);          %first session
                cat_Resp = squeeze(obj.catTrials(Resp));            % [total frames x cells]
                if strcmp(stim, 'PDG')
                    exidx = logical(repmat([ones(1, 40) zeros(1, 20)], 1, 96));
                    cat_Resp(exidx, :) = [];                        %excise offtime for PDG
                end
                zScore = mean(cat_Resp, 1)./std(cat_Resp, [], 1);

            case 'average'
                switch stim
                case 'NatMov'
                    Resp = obj.RespData.(stim).RespMat_onTime;          %average across sessions
                case 'PDG'
                    Resp = obj.RespData.(stim).RespMat_Full;
                end

                if strcmp(test, 'reliability')
                    reliability = nanmean(obj.StabilityData.(stim).CCs, 1);     %average across sessions
                end

                for kk = 1:obj.num_sessions                         %average across sessions
                    curr_Resp = permute(Resp(:, :, :, kk), [2 3 4 1]);
                    cat_Resp = squeeze(obj.catTrials(curr_Resp));            % [total frames x cells]
                    if strcmp(stim, 'PDG')
                        exidx = logical(repmat([ones(1, 40) zeros(1, 20)], 1, 96));
                        cat_Resp(exidx, :) = [];                        %excise offtime for PDG
                    end
                    zScore(:, kk) = nanmean(cat_Resp, 1)./std(cat_Resp, [], 1);
                end
                zScore = nanmean(zScore, 2);
            end
            mean_RDI = nanmean(obj.StabilityData.(stim).RDI(:, good_cells), 1);

            switch test
                case 'responsivity'
                    testmetric = zScore;
                case 'reliability'
                    testmetric = reliability;
            end

            manset = questdlg('Force axes limits?', 'Query', 'Yes', 'No', 'No');
            if strcmp(manset, 'Yes')
                xlimits = input('Enter x-axis limits\n');
                ylimits = input('Enter y-axis limits\n');
            end

            edges = prctile(testmetric(good_cells), 0:binsize:100);
            binned_zScore = discretize(testmetric(good_cells), edges);
            for ee = 1:length(edges)-1
                binned_RDI(ee) = mean(mean_RDI(binned_zScore == ee));
                binned_RDI_SEM(ee) = std(mean_RDI(binned_zScore == ee))/sqrt(sum(binned_zScore == ee));
                midpoints(ee) = (edges(ee)+edges(ee+1))/2;
            end

            mdl = fitlm(midpoints', binned_RDI');            % get confidence intervals for plotting later
            Xnew = linspace(min(midpoints)-0.5, max(midpoints)+0.5, 1000)';
            [~, CI] = predict(mdl, Xnew);

            [f, S] = polyfit(midpoints, binned_RDI, 1);
            rsq = 1 - (S.normr/norm(binned_RDI - mean(binned_RDI)))^2;

            %perform ranksum test on 1st and 4th quartiles
            top_RDI = mean_RDI(testmetric(good_cells) >= prctile(testmetric(good_cells), 75));
            bottom_RDI = mean_RDI(testmetric(good_cells) <= prctile(testmetric(good_cells), 25));
            sample_size = length(top_RDI);
            iterations = 1000;
            ps = [];
            for tt = 1:iterations
                bootsample = datasample(bottom_RDI, sample_size);
                ps(tt) = ranksum(top_RDI, bootsample);
            end
            p_boot = mean(ps);
            % fprintf('bootstrapped ranksum p = %.2e\n', p_boot);

            [p_thresh, ~, stats] = ranksum(top_RDI, bottom_RDI);

            % move outliers to boundaries of the plot if desired
            if exist('xlimits')
                if ~isempty(xlimits) 
                    testmetric(testmetric < xlimits(1)) = xlimits(1);
                    testmetric(testmetric > xlimits(2)) = xlimits(2);
                    midpoints(midpoints < xlimits(1)) = xlimits(1);
                    midpoints(midpoints > xlimits(2)) = xlimits(2);
                end
            end
            if exist('ylimits')
                if ~isempty(ylimits)
                    mean_RDI(mean_RDI < ylimits(1)) = ylimits(1);
                    mean_RDI(mean_RDI > ylimits(2)) = ylimits(2);
                    binned_RDI(binned_RDI < ylimits(1)) = ylimits(1);
                    binned_RDI(binned_RDI > ylimits(2)) = ylimits(2);
                    binned_RDI_SEM(binned_RDI_SEM < ylimits(1)) = ylimits(1);
                    binned_RDI_SEM(binned_RDI_SEM > ylimits(2)) = ylimits(2);
                end
            end
            figure
            scatter(testmetric(good_cells), mean_RDI, 'filled', 'MarkerFaceAlpha', 0.4);
            hold on
            errorbar(midpoints, binned_RDI, binned_RDI_SEM, 'vertical', 'LineStyle', 'none');
            scatter(midpoints, binned_RDI, 'filled');
            refline(f(1), f(2))
            plot(Xnew, CI(:, 1));
            plot(Xnew, CI(:, 2));
            [r_binned, p_binned] = corr(midpoints', binned_RDI', 'Type', 'Spearman');
            title(sprintf('r = %.2f, p = %.2e, rsq = %.2f', r_binned, p_binned, rsq))
            xlabel('testmetric');
            ylabel('Average RDI');
            if logscale == 1
                set(gca, 'xscale', 'log')
            end
            axis square
            if exist('xlimits')
                if ~isempty(xlimits) 
                    xlim([xlimits(1) xlimits(2)])
                end
            else
                xlim([min(testmetric(good_cells)) max(testmetric(good_cells))])
            end
            if exist('ylimits')
                if ~isempty(ylimits)
                    ylim([ylimits(1) ylimits(2)])
                end
            end

            figure
            hold on
            errorbar(midpoints, binned_RDI, binned_RDI_SEM, 'vertical', 'LineStyle', 'none');
            scatter(midpoints, binned_RDI, 'filled');
            f = polyfit(midpoints, binned_RDI, 1);
            % xlim([-0.2 1])
            
            mdl = fitlm(midpoints', binned_RDI');
            Xnew = linspace(min(midpoints)-0.2, max(midpoints)+0.2, 1000)';
            [~, CI] = predict(mdl, Xnew);
            plot(Xnew, CI(:, 1));
            plot(Xnew, CI(:, 2));
            xlabel('testmetric (binned)');
            ylabel('Average RDI (binned)')
            refline(f(1), f(2));
            xlim([min(midpoints)-0.2 max(midpoints)+0.2])
            axis square

            
            
            figure
            boxplot(cat(2, bottom_RDI, top_RDI), [zeros(1, length(bottom_RDI)) ones(1, length(top_RDI))]);
            title(sprintf('p = %.2e', p_thresh))
            xlabel('testmetric');
            ylabel('Average RDI');
            set(gcf, 'Position', [400 400 200 500]);

            % figure
            % hold on
            % histogram(top_RDI, 10, 'FaceColor', 'g', 'FaceAlpha', 0.3);
            % histogram(bottom_RDI, 10, 'FaceColor', 'r', 'FaceAlpha', 0.1);

            % figure 
            % scatter(zScore(good_cells), reliability(good_cells), 'filled');
            % axis square
            % [r, p] = corr(zScore(good_cells)', reliability(good_cells)', 'Type', 'Spearman');
            % title(sprintf('r = %.2f, p = %.2e', r, p))
            % xlabel('responsivity')
            % ylabel('reliability')

            plotdata.scatter.x = testmetric(good_cells);
            plotdata.scatter.y = mean_RDI;
            plotdata.scatter.xbinlocs = midpoints;
            plotdata.scatter.ybinmeans = binned_RDI;
            plotdata.scatter.shadedx = Xnew;
            plotdata.scatter.shadedy = CI(:, 2);
            plotdata.scatter.lineslope = f(1);
            plotdata.scatter.lineyint = f(2);
            plotdata.box.firstquartile = bottom_RDI;
            plotdata.box.lastquartile = top_RDI;
        end

        function waveformdata = extractWaveforms(obj, stim, min_eventsize, pvalue, min_sessionthresh)
            % waveformdata output is a struct with 3 fields:    1) event_vector: [cells x frames] logical vectors denoting the event periods across the recording for every cell 
            %                                                   2) waveform_ts:  [1 x cells] cell array of [num_waves x 2] arrays containing start (first column) and end (second column) timestamps for every event for every cell
            %                                                   3) waveforms:    [1 x cells] cell array of [1 x num_waves] cell arrays of [repeats x frames x sessions] vectors of extracted waveforms of every event for every cell

            waveformdata = struct();            % initialize output struct
            waveformdata.event_vectors = [];
            waveformdata.waveform_ts = {};
            waveformdata.waveforms = {};

            Resp = obj.RespData.(stim).RespMat_Full;

            Resp_sp = obj.spikeResp.(stim).RespMat_Full;

            switch stim
                case 'NatMov'
                    F0_idx = 31:50;     % offtime frames for determining F0
                    off_idx = 1:50;
                case 'PDG'
                    F0_idx = logical(repmat([zeros(1, 20) ones(1, 20) zeros(1, 20)], 1, 12));
                    off_idx = logical(repmat([ones(1, 40) zeros(1, 20)], 1, 12));
            end

            offResp = squeeze((mean(Resp(:, F0_idx, :, :), 2)));           %average offtime response on every trial
            offResp_sp = squeeze((mean(Resp_sp(:, F0_idx, :, :), 2)));         % for spike data
            
            cellnums = obj.num_cells_bymouse;         % length of this vector is how many fields there are, element values are how many cells per field (ALL NEURONS)
            if isempty(cellnums)                % in case we're only doing one field
                cellnums = obj.num_cells;
            end
            
            sesnums = obj.num_sessions_bymouse;
            if isempty(sesnums)                 % in case we're only doing one field
                sesnums = obj.num_sessions;
            end
            num_fields = length(cellnums);

            % need to perform event detection and waveform extraction on each field separately, and add them to the growing fields of the output struct
            % thus, need to run the following in blocks of cells corresponding to cells in each field
            running_sum = 0;        % tracks number of neurons we've been through so far
            for zz = 1:num_fields
                fprintf('Extracting waveforms for field %d out of %d\n', zz, num_fields);
                curr_cellidcs = running_sum+1:running_sum + cellnums(zz);           % cell indices for current field
                running_sum = running_sum + cellnums(zz);

                curr_sesidcs = 1:sesnums(zz);                                         % session indicies for current field

                currfield_Resp = Resp(:, :, curr_cellidcs, curr_sesidcs);
                currfield_Resp_sp = Resp_sp(:, :, curr_cellidcs, curr_sesidcs);
                currfield_offResp = offResp(:, curr_cellidcs, curr_sesidcs);
                currfield_offResp_sp = offResp_sp(:, curr_cellidcs, curr_sesidcs);

                event_vectors = eventDetector(currfield_Resp, currfield_offResp, currfield_Resp_sp, currfield_offResp_sp, min_eventsize, pvalue, min_sessionthresh);
                event_vectors(:, off_idx) = 0;                                          % discard any events found during the offtime
                waveformdata.event_vectors = cat(1, waveformdata.event_vectors, event_vectors);     % add to output

                waveform_ts = cell(1, cellnums(zz));
                waveforms = cell(1, cellnums(zz));
                for ii = 1:cellnums(zz)            % iterate through neurons of this field, use event vector to extract remaining output fields
                    %some cleanup to deal with residual event frames resulting from discarding offtime event frames (this is the same function that occurs in the original event detection function)
                    curr_vector = event_vectors(ii, :);
                    switches = diff(curr_vector);               
                    switch_idcs = find(switches == 1);          % indices (respective to 'switches') where curr_vector goes from 0 to 1
                    for qq = 1:length(switch_idcs)              % iterate through every index where this occurs, count the number of contiguous 0s (which corresponds to no element value change, in this case we're counting 1s)
                        curr_idx = switch_idcs(qq);
                        ele_counter = 1;
                        end_flag = 0;
                        while ~end_flag
                            if (curr_idx + ele_counter) <= length(switches)             % check in case we reach the end of the vector
                                curr_val = switches(curr_idx + ele_counter);

                                if curr_val == 0                            % if we encounter a zero, increment our contiguity counter
                                    ele_counter = ele_counter + 1;
                                else
                                    end_flag = 1;                           % if we encounter not a zero, it means we've gone from 1 to 0 in curr_vector and the string of contiguous 1s is over, break the while loop
                                end
                            else
                                end_flag = 1;
                            end
                        end
                        if ele_counter < min_eventsize                         % number of times we incremented our contiguity counter corresponds to the length of the 'event' we're looking at
                            curr_vector(curr_idx+1:curr_idx+ele_counter) = 0;       % if it's less than a minimum length, set all the values of this 'event' to 0
                        end
                    end
                    event_vectors(ii, :) = curr_vector;

                    %find event timestamps
                    startstops = diff(event_vectors(ii, :));          % 1s (+1) indicate where events start, -1s indicate where events end

                    starts = find(startstops == 1)+1;
                    stops = find(startstops == -1);

                    % some cleanup to account for events that start or stop on the first or last frames
                    if length(starts) > length(stops)       % case in which event starts on first frame
                        stops(end+1) = size(Resp, 2);
                    elseif length(stops) > length(starts)       % case in which event stops on last frame
                        starts = [1 starts];
                    elseif length(starts) > 1                   % case in which both of the above happen in the same neuron
                        if starts(1) > stops(1)
                            stops(end+1) = size(Resp, 2);
                            starts = [1 starts];
                        end
                    end

                    %some cleanup to account for single rogue event frames (because eventDetector still has bugs)
                    ww = 1;
                    while ww <= length(starts)
                        if stops(ww) - starts(ww) == 0
                            event_vectors(ii, starts(ww)) = 0;
                            starts(ww) = [];
                            stops(ww) = [];
                        else
                            ww = ww + 1;
                        end
                    end

                    waveform_ts{ii} = cat(2, starts', stops');

                    for ww = 1:length(starts)                       % extract every waveform for this cell
                        waveforms{ii}{ww} = squeeze(currfield_Resp(:, starts(ww):stops(ww), ii, :));                        % currently outputting original data, not normalized
                    end
                end
                waveformdata.waveform_ts = cat(2, waveformdata.waveform_ts, waveform_ts);
                waveformdata.waveforms = cat(2, waveformdata.waveforms, waveforms);
            end
        end

        
        % Submethods
        function color = getColorMap(obj, stim_type)
            switch stim_type
                case 'PDG'
                    color = [0.0 0.7 0.7];
                case 'NatMov'
                    color = [0.8 0.5 0.5];
                case 'MultiMov'
                    call_count = obj.MultiMov_called;
                    color = [0.7+0.1*(call_count-1) 0.4+0.2*(call_count-1) 0.4+0.1*(call_count-1)];
                    %some smart shit
                case 'PDGcont'
                    color = [0.4 0.6 0.9];
                case 'NatMovdisc'
                    color = [0.9 0.3 0.3];
            end
        end

        function [RDIdata, n] = extractRDI(obj, stim_type)
            use_cells_stretched = repmat(obj.getUse_cells', 1, obj.num_sessions-1);
            use_sessions = obj.use_sessions_full(:, 2:obj.num_sessions);
            use_sessions(isnan(use_sessions)) = 0;
            combined_indicator = use_cells_stretched & use_sessions;          %combined the cell indicator and the session indicator for easy indexing below
            n = sum(combined_indicator, 1);
            
            Stabdata = obj.getStabilityData;
            this_stim_stabdata = Stabdata.(stim_type);

            if ~iscell(this_stim_stabdata.RDI)
                RDI = this_stim_stabdata.RDI;
                RDI(RDI > 1) = 1;
                RDI(RDI < -1) = -1;
                RDIdata.RDI = RDI;
                RDIdata.RDI_included_scatter = RDI;
                RDIdata.RDI_included_scatter(~combined_indicator') = NaN;
                for jj = 1:obj.num_sessions-1
                    RDIdata.RDI_avg(jj) = nanmean(RDI(jj, combined_indicator(:, jj)));
                    RDIdata.RDI_SEM(jj) = nanstd(RDI(jj, combined_indicator(:, jj)), [], 2)/sqrt(sum(combined_indicator(:, jj)));
                end
            else
                for kk = 1:length(this_stim_stabdata.RDI)
                    RDI = this_stim_stabdata.RDI{kk};
                    RDI(RDI > 1) = 1;
                    RDI(RDI < -1) = -1;
                    RDIdata.RDI{kk} = RDI;
                    RDIdata.RDI_included_scatter{kk} = RDI;
                    RDIdata.RDI_included_scatter{kk}(~combined_indicator') = NaN;
                    for jj = 1:obj.num_sessions-1
                        RDIdata.RDI_avg{kk}(jj) = nanmean(RDI(jj, combined_indicator(:, jj)));
                        RDIdata.RDI_SEM{kk}(jj) = nanstd(RDI(jj, combined_indicator(:, jj)), [], 2)/sqrt(sum(combined_indicator(:, jj)));
                    end
                end
            end
        end

        % Getters
        function RDIdata = getRDIdata(obj)
            RDIdata = obj.RDIdata;
        end
        
        function setRDIdata(obj, input)
            obj.RDIdata = input;
        end

        function setSpikedata(obj, input)
            %input is a StabilityAnalyzer object with spike RespData
            obj.spikeResp = input.RespData;
        end
    end
end