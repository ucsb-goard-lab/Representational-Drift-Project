classdef EnsembleAnalyzer < Analyzer
    
    properties
        Spontdata
        spont_ref
    end
    
    
    methods
        
        function obj = EnsembleAnalyzer(varargin)
            %Initilize the object with desired data
            %Inputs can be any of the above properties, in any order, as long as their names exactly match 
            %If no inputs are provided, you can select them from your directory, also in any order
            proplist =  properties(obj);
            if nargin == 0
                disp('Select all data files')
                dataNames = uigetfile('.mat', 'MultiSelect', 'on');
                for zz = 1:length(dataNames)
                    subroutine_progressbar(zz/length(dataNames));
                    [~, loc] = ismember(erase(dataNames{zz}, '.mat'), proplist);
                    obj.(proplist{loc}) = importdata(dataNames{zz});
                end
                close all
            else
                for ii = 1:nargin
                    varname = inputname(ii);
                    [~, loc] = ismember(varname, proplist);
                    obj.(proplist{loc}) = varargin{ii};
                end
            end 
            if ~isempty(obj.RespData)
                try
                    obj.num_sessions = size(obj.RespData.PDG.RespMat_Full, 4);
                    obj.num_cells = size(obj.RespData.PDG.RespMat_Full, 3);
                catch
                    try
                        obj.num_sessions = size(obj.RespData.NatMov.RespMat_Full, 4);        %add more try catches here for additional stimuli
                        obj.num_cells = size(obj.RespData.NatMov.RespMat_Full, 3);
                    catch
                        obj.num_sessions = size(obj.RespData.Scenes.RespMat_Full, 4);
                        obj.num_cells = size(obj.RespData.Scenes.RespMat_Full, 3);
                    end
                end
            end
        end
                
        function importSpontaneous(obj)
            fprintf('Select all spontaneous activity datafiles...\n');
            filenames = uigetfile('.mat', 'MultiSelect', 'on');
            if iscell(filenames)
                for ii = 1:length(filenames)
                    currdata = importdata(filenames{ii});
                    RespMat(1, :, :, ii) = transpose(currdata.DFF);
                end
            else
                currdata = importdata(filenames);
                RespMat(1, :, :, 1) = transpose(currdata.DFF);
            end
            
            obj.RespData.Spont = structPacker([], RespMat, 'RespMat_Full');
            if isempty(obj.num_sessions)
                obj.num_sessions = size(RespMat, 4);
            end
            if isempty(obj.num_cells)
                obj.num_cells = size(RespMat, 3);
            end
        end
        
        function setSpontReference(obj, ref)
            obj.spont_ref = ref;
        end

        
        function corr_mat = SignalCorr(obj, RespMat1, RespMat2)          % CC between trial averaged activity 
            num_cells = size(RespMat1, 3);
            corr_mat = zeros(num_cells, num_cells);
            if nargin < 3                                       %if only one resp matrix provided, use all trials for both cells
                avg_traces = squeeze(mean(RespMat1, 1));
                for ii = 1:num_cells
                    disp(ii);
                    for jj = 1:num_cells
                        corr_mat(ii, jj) = corr(avg_traces(:, ii), avg_traces(:, jj));
                    end
                end
            else                                                % if multiple matrices provided, correlate using both where iterating ii = first respmat and iterating jj = second respmat
                avg_traces1 = squeeze(mean(RespMat1, 1));
                avg_traces2 = squeeze(mean(RespMat2, 1));
                for ii = 1:num_cells
                    disp(ii);
                    for jj = 1:num_cells
                        corr_mat(ii, jj) = corr(avg_traces1(:, ii), avg_traces2(:, jj));
                    end
                end
            end
        end

        function corr_mat = WholeCorr(obj, RespMat, num_frames)          % CC between raw full traces
            num_sessions = obj.num_sessions;

            qual = obj.getUse_cells;
            sessions = obj.getUse_sessions;
            sessions(isnan(sessions)) = 1;
            present = sum(sessions, 2) == num_sessions;
            good_cells = qual' & present;

            RespMat = RespMat(:, :, good_cells);

            traces_cat = [];
            for rr = 1:size(RespMat, 1)
                traces_cat = cat(1, traces_cat, squeeze(RespMat(rr, :, :))); 
            end
            
            corr_mat = zeros(sum(good_cells), sum(good_cells));
            for ii = 1:size(RespMat, 3)
                for jj = 1:size(RespMat, 3)
                    corr_mat(ii, jj) = corr(traces_cat(1:num_frames, ii), traces_cat(1:num_frames, jj));
                end
            end
        end

        function corr_mat = NoiseCorr(obj, RespMat, num_frames, mean_activity)          % CC between de-meaned full traces    
            num_sessions = obj.num_sessions;

            qual = obj.getUse_cells;
            sessions = obj.getUse_sessions;
            sessions(isnan(sessions)) = 1;
            present = sum(sessions, 2) == num_sessions;
            good_cells = qual' & present;

            RespMat = RespMat(:, :, good_cells);

            if nargin == 4
                avg_traces = mean_activity(:, good_cells);
            else
                avg_traces = squeeze(mean(RespMat, 1));
            end

            for rr = 1:size(RespMat, 1)
                noise_traces(rr, :, :) = squeeze(RespMat(rr, :, :)) - avg_traces;
            end

            noise_traces_cat = [];
            for rr = 1:size(noise_traces, 1)
                noise_traces_cat = cat(1, noise_traces_cat, squeeze(noise_traces(rr, :, :))); 
            end
            
            corr_mat = zeros(sum(good_cells), sum(good_cells));
            for ii = 1:size(RespMat, 3)
                for jj = 1:size(RespMat, 3)
                    corr_mat(ii, jj) = corr(noise_traces_cat(1:num_frames, ii), noise_traces_cat(1:num_frames, jj));
                end
            end
        end

        function corr_mat = NoiseCorr_v2(obj, RespMat, num_trials, mean_activity)     % trials are de-meaned, CC found for every trial, averaged across trials   
            num_sessions = obj.num_sessions;

            qual = obj.getUse_cells;
            sessions = obj.getUse_sessions;
            sessions(isnan(sessions)) = 1;
            present = sum(sessions, 2) == num_sessions;
            good_cells = qual' & present;

            RespMat = RespMat(:, :, good_cells);

            if nargin == 3
                avg_traces = mean_activity(:, good_cells);
            else
                avg_traces = squeeze(mean(RespMat, 1));
            end

            for rr = 1:size(RespMat, 1)
                noise_traces(rr, :, :) = squeeze(RespMat(rr, :, :)) - avg_traces;
            end
            
            corr_mat = zeros(sum(good_cells), sum(good_cells));
            for ii = 1:size(RespMat, 3)
                for jj = 1:size(RespMat, 3)
                    temp = zeros(1, num_trials);
                    for rr = 1:num_trials
                        temp(rr) = corr(squeeze(noise_traces(rr, :, ii))', squeeze(noise_traces(rr, :, jj))');
                    end
                    corr_mat(ii, jj) = mean(temp);
                end
            end
        end

        function binned_Resp = binner(obj, Resp, stim)
            bin_size = 20;          %frames
            interval_size = 40;      %frames in between
            chunksize = bin_size + interval_size;

            num_reps = size(Resp, 1);
            num_frames = size(Resp, 2);
            num_cells = size(Resp, 3);

            %if natmov detected, restructure first
            if num_reps == 30 || num_reps == 20
                restruct_resp = zeros(num_reps/2, num_frames*2, num_cells);
                for rr = 1:num_reps/2
                    newtrial = cat(2, Resp(rr*2-1, :, :), Resp(rr*2, :, :));
                    restruct_resp(rr, :, :) = newtrial;
                end
                Resp = restruct_resp;
            end

            num_frames = size(Resp, 2);
            num_bins = ceil(num_frames/chunksize);

            binned_Resp = zeros(size(Resp, 1), num_bins, num_cells);
            for ck = 1:num_bins
                if strcmp(stim, 'NatMov')
                    curr_idx = (ck-1)*chunksize+1:(ck-1)*chunksize+bin_size;
                elseif strcmp(stim, 'PDG')
                    curr_idx = ck*chunksize-20+1:ck*chunksize;
                end
                binned_Resp(:, ck, :) = mean(Resp(:, curr_idx, :), 2);
            end
        end

        function corr_mat = NoiseCorr_Montijn(obj, Resp)                % activity is binned, CC found between trial-to-trial vectors of binned activity, averaged across bins

            num_sessions = obj.num_sessions;

            qual = obj.getUse_cells;
            sessions = obj.getUse_sessions;
            sessions(isnan(sessions)) = 1;
            present = sum(sessions, 2) == num_sessions;
            good_cells = qual' & present;

            Resp = Resp(:, :, good_cells);          % do this here to save time instead of doing the whole matrix 

            num_cells = size(Resp, 3);
            corr_mat = zeros(num_cells, num_cells);
            num_bins = size(Resp, 2);

            for ii = 1:num_cells
                % disp(ii);
                for jj = 1:num_cells
                    curr_corr = zeros(1, num_bins);
                    for bb = 1:num_bins
                        curr_corr(bb) = corr(Resp(:, bb, ii), Resp(:, bb, jj));         % find correlation of trial-to-trial response vector for every bin
                    end
                    corr_mat(ii, jj) = mean(curr_corr);             %mean across bins
                end
            end
        end

        function dist = popnavgDistribution(obj, RespMat, session)
            iterations = 100;
            qual = obj.getUse_cells;
            present = obj.getUse_sessions;
            good_cells = qual' & present;           % [cells x sessions] good and present on session kk
            numReps = size(RespMat, 1);
            cat_RespMat = [];

            if numReps == 1
                cat_RespMat = squeeze(RespMat);
            else
                for rr = 1:numReps
                    cat_RespMat = cat(1, cat_RespMat, squeeze(RespMat(rr, :, :)));           %concatenate all reps for full recording
                end
            end

            dist = zeros(iterations, size(cat_RespMat, 1));

            for hh = 1:iterations
                hh
                choosefrom = find(good_cells(:, session));      % use all good cells present on this session
                shuffled = randperm(length(choosefrom));                % shuffled list of indices of choosefrom
                idx_list = shuffled(1:floor(length(choosefrom)/2));          % random half of neurons

                dist(hh, :) = nanmean(cat_RespMat(:, idx_list), 2);
            end
        end

        function out = getPDGofftimecat(obj)
            in = obj.RespData.PDG.RespMat_Full;
            offtime = 40;       %frames
            ontime = 20;
            num_ori = 12;
            frame_idx = [ones(1, offtime) zeros(1, ontime)];
            frame_idx = logical(repmat(frame_idx, 1, num_ori));
            out = in(:, frame_idx, :, :);
        end
            
    end
end