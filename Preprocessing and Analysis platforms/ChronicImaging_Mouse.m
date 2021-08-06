classdef ChronicImaging_Mouse < General_Processor
    
    properties
        OrientationData     % data structure containing orientation tuning data
        RespData            % data structure containing neural response matrices for all stimuli
        StabilityData       % data structure containing reliability and stability data for all stimuli
        RFdata              % data structure containig receptive field data 
        RoiINFO             % data structure containing visual responsiveness and ROI quality and presence data
        Vis_cellInspection_maps     % saved cell map data for visual inspection software
        PERCENTILE = 99.9;  % significance threshold for visual responsiveness procedure    
        subsample_flag      % whether or not to subsample trials for data in use 
    end

    properties (Access = private)
        PDG_rel_shuffletest double
        NatMov_rel_shuffletest double
        ZSCORE = 1.5;
        ALPHA = 0.01; 
        subsample           % value for subsampling trials
    end
    
    methods
        
        function obj = ChronicImaging_Mouse()
            obj.OrientationData = struct();
            obj.RespData = struct();
            obj.StabilityData = struct();
            obj.RFdata = struct();
            obj.RoiINFO = struct();
            obj.Vis_cellInspection_maps = struct();
        end
        
        function importObjects(obj, varargin)
            for ii = 1:2:length(varargin)
                curr_import = varargin{ii};
                curr_dataobj = varargin{ii+1};
                currRespData = curr_dataobj.processor.getRespData;
                if ~strcmp(curr_import, 'RF')
                    currStabilityData = curr_dataobj.processor.getStabilityData;
                    obj.RespData.(curr_import) = structPacker([], currRespData.RespMat_Full, 'RespMat_Full');
                    try
                        obj.RespData.(curr_import) = structPacker(obj.RespData.(curr_import), currRespData.RespMat_onTime, 'RespMat_onTime');
                    catch
                    end
                    try
                        obj.RespData.(curr_import) = structPacker(obj.RespData.(curr_import), currRespData.SceneResp_norm, 'SceneResp_norm');
                    catch
                    end
                    obj.StabilityData.(curr_import) = structPacker([], currStabilityData.CCs, 'CCs', currStabilityData.RDI, 'RDI',...
                        currStabilityData.CC_ws, 'CC_ws', currStabilityData.CC_bs, 'CC_bs', currStabilityData.RDI_control, 'RDI_control');
                    try
                        OrientationData = curr_dataobj.processor.getOrientationData;
                        if ~isempty(fieldnames(OrientationData))
                            obj.OrientationData = curr_dataobj.processor.getOrientationData;
                        end
                    catch
                    end     
                else
                    obj.RespData.RF = structPacker([], currRespData.altitude, 'altitude', ...
                        currRespData.azimuth, 'azimuth');
                    obj.RFdata.RespData = obj.RespData.RF;
                    obj.RFdata.RFmapping_results = curr_dataobj.processor.getRFdata;
                end
            end
            if ~isempty(obj.RespData) && isempty(fieldnames(obj.RFdata))
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
            else
                obj.num_sessions = size(obj.RFdata.RFmapping_results.isSpatiallyTuned, 1);
                obj.num_cells = size(obj.RFdata.RFmapping_results.isSpatiallyTuned, 2);
             %for importing preprocessor objects
            end
        end

        function importData(obj, weeks, varargin)
            %Initilize the object with desired data
            %Inputs 2:n can be any of the above properties, in any order, as long as their names exactly match 
            %If no inputs are provided, you can select them from your directory, also in any order
            proplist =  properties(obj);
            if nargin < 3
                disp('Select all data files (RespData, RoiINFO, & StabilityData)')
                dataNames = uigetfile('.mat', 'MultiSelect', 'on');
                if iscell(dataNames)
                    for zz = 1:length(dataNames)
                        % subroutine_progressbar(zz/length(dataNames));
                        [~, loc] = ismember(erase(dataNames{zz}, '.mat'), proplist);
                        obj.(proplist{loc}) = importdata(dataNames{zz});
                    end
                else
                    [~, loc] = ismember(erase(dataNames, '.mat'), proplist);            % in case only loading one file
                    obj.(proplist{loc}) = importdata(dataNames);
                end
            else
                for ii = 1:nargin
                    varname = inputname(ii);
                    [~, loc] = ismember(varname, proplist);
                    obj.(proplist{loc}) = varargin{ii};
                end
            end 

            if nargin < 2
                weeks = [];
            end

            if ~isempty(weeks)
                aligner = Analyzer();
                aligner.setWeeks(weeks);

                if length(fieldnames(obj.RespData)) ~= 0
                    obj.RespData.PDG.RespMat_Full = aligner.alignWeeks(obj.RespData.PDG.RespMat_Full, 4, 3, 'RespData');
                    try
                    obj.RespData.NatMov.RespMat_Full = aligner.alignWeeks(obj.RespData.NatMov.RespMat_Full, 4, 3, 'RespData');
                    obj.RespData.NatMov.RespMat_onTime = aligner.alignWeeks(obj.RespData.NatMov.RespMat_onTime, 4, 3, 'RespData');
                    catch
                    end
                end

                if length(fieldnames(obj.StabilityData)) ~= 0
                    obj.StabilityData.PDG.RDI = aligner.alignWeeks(obj.StabilityData.PDG.RDI, 1, 2, 'StabilityData', -1);           % correction factor of -1 due to first index being week 2
                    obj.StabilityData.PDG.CCs = aligner.alignWeeks(obj.StabilityData.PDG.CCs, 1, 2, 'StabilityData');
                    obj.StabilityData.PDG.CC_ws = aligner.alignWeeks(obj.StabilityData.PDG.CC_ws, 1, 2, 'StabilityData', -1);
                    obj.StabilityData.PDG.CC_bs = aligner.alignWeeks(obj.StabilityData.PDG.CC_bs, 1, 2, 'StabilityData', -1);
                    try
                    obj.StabilityData.NatMov.RDI = aligner.alignWeeks(obj.StabilityData.NatMov.RDI, 1, 2, 'StabilityData', -1);           % correction factor of -1 due to first index being week 2
                    obj.StabilityData.NatMov.CCs = aligner.alignWeeks(obj.StabilityData.NatMov.CCs, 1, 2, 'StabilityData');
                    obj.StabilityData.NatMov.CC_ws = aligner.alignWeeks(obj.StabilityData.NatMov.CC_ws, 1, 2, 'StabilityData', -1);
                    obj.StabilityData.NatMov.CC_bs = aligner.alignWeeks(obj.StabilityData.NatMov.CC_bs, 1, 2, 'StabilityData', -1);                        
                    catch
                    end
                end

                if length(fieldnames(obj.OrientationData)) ~= 0
                    obj.OrientationData.isTuned = aligner.alignWeeks(obj.OrientationData.isTuned, 1, 2, 'OrientationData');     
                    obj.OrientationData.osiMat = aligner.alignWeeks(obj.OrientationData.osiMat, 1, 2, 'OrientationData'); 
                    obj.OrientationData.zScoreMat = aligner.alignWeeks(obj.OrientationData.zScoreMat, 1, 2, 'OrientationData');
                    obj.OrientationData.oriResp = aligner.alignWeeks(obj.OrientationData.oriResp, 1, 2, 'OrientationData');
                    obj.OrientationData.oriPref = aligner.alignWeeks(obj.OrientationData.oriPref, 1, 2, 'OrientationData');
                    obj.OrientationData.osiAvg = aligner.alignWeeks(obj.OrientationData.osiAvg, 2, 0, 'OrientationData');
                    obj.OrientationData.zScoreAvg = aligner.alignWeeks(obj.OrientationData.zScoreAvg, 2, 0, 'OrientationData');
                    obj.OrientationData.osiSEM = aligner.alignWeeks(obj.OrientationData.osiSEM, 2, 0, 'OrientationData');
                    obj.OrientationData.zScoreSEM = aligner.alignWeeks(obj.OrientationData.zScoreSEM, 2, 0, 'OrientationData');   
                end

                if length(fieldnames(obj.RoiINFO)) ~= 0
                    obj.RoiINFO.presence = aligner.alignWeeks(obj.RoiINFO.presence, 2, 1, 'RoiINFO');     
                end

                if length(fieldnames(obj.RFdata)) ~= 0
                    obj.RFdata.RFmapping_results.alt_fit = aligner.alignWeeks(obj.RFdata.RFmapping_results.alt_fit, 1, 2, 'RFdata');     
                    obj.RFdata.RFmapping_results.azi_fit = aligner.alignWeeks(obj.RFdata.RFmapping_results.azi_fit, 1, 2, 'RFdata'); 
                    obj.RFdata.RFmapping_results.alt_pref = aligner.alignWeeks(obj.RFdata.RFmapping_results.alt_pref, 1, 2, 'RFdata');
                    obj.RFdata.RFmapping_results.azi_pref = aligner.alignWeeks(obj.RFdata.RFmapping_results.azi_pref, 1, 2, 'RFdata');
                    obj.RFdata.RFmapping_results.alt_p = aligner.alignWeeks(obj.RFdata.RFmapping_results.alt_p, 1, 2, 'RFdata');
                    obj.RFdata.RFmapping_results.azi_p = aligner.alignWeeks(obj.RFdata.RFmapping_results.azi_p, 1, 2, 'RFdata');
                    obj.RFdata.RFmapping_results.isSpatiallyTuned = aligner.alignWeeks(obj.RFdata.RFmapping_results.isSpatiallyTuned, 1, 2, 'RFdata');
                    obj.RFdata.RFmapping_results.roi_centroids = aligner.alignWeeks(obj.RFdata.RFmapping_results.roi_centroids, 3, 1, 'RFdata');
                    obj.RFdata.RFmapping_results.alt_corr = aligner.alignWeeks(obj.RFdata.RFmapping_results.alt_corr, 2, 0, 'RFdata');
                    obj.RFdata.RFmapping_results.azi_corr = aligner.alignWeeks(obj.RFdata.RFmapping_results.azi_corr, 2, 0, 'RFdata');
                    obj.RFdata.RFmapping_results.isVisuallyResponsive_alt = aligner.alignWeeks(obj.RFdata.RFmapping_results.isVisuallyResponsive_alt, 2, 1, 'RFdata');
                    obj.RFdata.RFmapping_results.isVisuallyResponsive_azi = aligner.alignWeeks(obj.RFdata.RFmapping_results.isVisuallyResponsive_azi, 2, 1, 'RFdata');
                    obj.RFdata.RFmapping_results.RespMat_alt = aligner.alignWeeks(obj.RFdata.RFmapping_results.RespMat_alt, 5, 1, 'RFdata');
                    obj.RFdata.RFmapping_results.RespMat_azi = aligner.alignWeeks(obj.RFdata.RFmapping_results.RespMat_azi, 5, 1, 'RFdata');
                end
            end
            if ~isempty(fieldnames(obj.RespData))

                try
                    obj.num_sessions = size(obj.RespData.PDG.RespMat_Full, 4);
                    obj.num_cells = size(obj.RespData.PDG.RespMat_Full, 3);
                catch
                    try
                        obj.num_sessions = size(obj.RespData.Scenes.RespMat_Full, 4);
                        obj.num_cells = size(obj.RespData.Scenes.RespMat_Full, 3);
                    catch
                        obj.num_sessions = size(obj.RespData.NatMov.RespMat_Full, 4);        %add more try catches here for additional stimuli
                        obj.num_cells = size(obj.RespData.NatMov.RespMat_Full, 3);
                    end
                end
            else
                obj.num_sessions = size(obj.RFdata.RFmapping_results.isSpatiallyTuned, 1);
                obj.num_cells = size(obj.RFdata.RFmapping_results.isSpatiallyTuned, 2);
                %for directly importing processed data structs
            end
        end
        
        function evaluateROIs(obj)
            %This method returns vectors that describe cells usable for stability analysis based on different criteria
            %Provide the method with any of the strings listed below in any order to return their respective vectors
            try
                obj.num_cells = size(obj.StabilityData.PDG.PDG_CCs, 2);
                obj.num_sessions = size(obj.StabilityData.PDG.PDG_CCs, 1);
            catch
            end
            possible_criteria = {'PDG_isTuned', 'PDG_Responsive', 'NatMov_Responsive', 'MultiMov_Responsive',...
                                'Visual_Inspection_A', 'Visual_Inspection_upload'};
            [selected_criteria, ~] = listdlg('ListString', possible_criteria);
            use_criteria = possible_criteria(selected_criteria);
            CC_thresh_PDG = prctile(obj.StabilityData.PDG.CCs(1, :), 25);          % must be above lowest quartile
            try
                CC_thresh_Nat = prctile(obj.StabilityData.NatMov.CCs(1, :), 25);      
            catch
            end
            tuned_thresh = obj.num_sessions;
            if any(strcmp(use_criteria, 'PDG_isTuned'))
                criterion = 'isTuned';
                PDG_isTuned =  sum(obj.OrientationData.isTuned, 1) >= tuned_thresh;
                fprintf('n = %d cells out of %d, using %s method \n', sum(PDG_isTuned), obj.num_cells, criterion)
                obj.RoiINFO.PDG_isTuned = PDG_isTuned;
            end
            if any(strcmp(use_criteria, 'PDG_Responsive'))
                % criterion = 'PDG biased CC';
                obj.setSubsampleFlag('No');

                [PDG_rel_shuffletest, threshold_distribution] = obj.trialShuffleTest(obj.RespData.PDG.RespMat_Full, obj.StabilityData.PDG.CCs, 'No');   
                obj.PDG_rel_shuffletest = PDG_rel_shuffletest;
%                 keyboard
                if obj.num_sessions > 2
                    PDG_statistical_thresh = ttest(obj.StabilityData.PDG.CCs, mean(threshold_distribution)) & mean(obj.StabilityData.PDG.CCs, 1) > mean(threshold_distribution);        %is mean ws CC statistically greater than threshold value
                else
                    PDG_statistical_thresh = mean(obj.StabilityData.PDG.CCs, 1) > (mean(threshold_distribution) + 3*std(threshold_distribution)) & mean(obj.StabilityData.PDG.CCs, 1) > CC_thresh_PDG;
                end
                obj.RoiINFO.PDG_Responsive_shuffle = PDG_rel_shuffletest;   %without threshold 
                obj.RoiINFO.PDG_Responsive_thresh = PDG_statistical_thresh;

            end     
            if any(strcmp(use_criteria, 'NatMov_Responsive'))
                % criterion = 'NatMov biased CC';
                obj.setSubsampleFlag('No');
                [NatMov_rel_shuffletest, threshold_distribution] = obj.trialShuffleTest(obj.RespData.NatMov.RespMat_onTime, obj.StabilityData.NatMov.CCs, 'No');
                obj.NatMov_rel_shuffletest = NatMov_rel_shuffletest;
%                 keyboard
                if obj.num_sessions > 2
                    NatMov_statistical_thresh = ttest(obj.StabilityData.NatMov.CCs, mean(threshold_distribution)) & mean(obj.StabilityData.NatMov.CCs, 1) > mean(threshold_distribution);  %Cells responsive to NatMov
                else 
                    NatMov_statistical_thresh = mean(obj.StabilityData.NatMov.CCs, 1) > (mean(threshold_distribution) + 3*std(threshold_distribution)) & mean(obj.StabilityData.NatMov.CCs, 1) > CC_thresh_Nat;
                end
                obj.RoiINFO.NatMov_Responsive_shuffle = NatMov_rel_shuffletest;
                obj.RoiINFO.NatMov_Responsive_thresh = NatMov_statistical_thresh;
      
            end
            if any(strcmp(use_criteria, 'MultiMov_Responsive'))
                criterion = 'MultiMov biased CC';
                [NatMov_rel_shuffletest, threshold_distribution] = obj.trialShuffleTest(obj.RespData.MultiMov.RespMat_onTime{1}, obj.StabilityData.MultiMov.CCs{1});
                obj.NatMov_rel_shuffletest = NatMov_rel_shuffletest;
                NatMov_statistical_thresh = ttest(obj.StabilityData.MultiMov.CCs{1}, mean(threshold_distribution)) & mean(obj.StabilityData.MultiMov.CCs{1}, 1) > mean(threshold_distribution);  %Cells responsive to NatMov

                obj.RoiINFO.MultiMov_Responsive_shuffle = NatMov_rel_shuffletest;
                obj.RoiINFO.MultiMov_Responsive_thresh = NatMov_statistical_thresh;
            end
            if any(strcmp(use_criteria, 'Visual_Inspection_A'))
                savedata_flag = 0;
                preprocess_flag = 0;
                loaddata = questdlg('Load saved data?', 'Import data', 'Yes', 'No', 'No');
                if strcmp(loaddata, 'Yes')
                    savedata_flag = 1;
                    cellinfofile = uigetfile('.mat');
                    cellinfo = importdata(cellinfofile);
                    startcell = cellinfo.lastCell;
                end
                if isempty(fieldnames(obj.Vis_cellInspection_maps))
                    preprocess_flag = 1;
                end
                if preprocess_flag == 0 
                    overwrite_flag = questdlg('Overwrite preprocessed data?', 'Overwrite', 'Yes', 'No', 'Yes');
                    if strcmp(overwrite_flag, 'Yes')
                        preprocess_flag = 1;
                    end
                end
                if preprocess_flag == 1
                    fprintf('Select NatMov maps file.\n');
                    filename = uigetfile('.mat');
                    natmov_maps = importdata(filename);

                    fprintf('Select PDG maps file.\n');
                    filename = uigetfile('.mat');
                    pdg_maps = importdata(filename);

                    fprintf('Select reference data file for cell masks.\n');
                    reference_filename = uigetfile('.mat');   
                    reference_data = importdata(reference_filename);
                    cellMasks = reference_data.cellMasks;

                    temp = inputdlg('Enter box size (pixels), default is 40 for 2x zoom.');
                    box_size = str2double(temp{1});
                    
                    fprintf('Processing NatMov maps.\n');
                    natmov_cellInspection = obj.visualCellInspection_preprocess(natmov_maps, box_size, cellMasks);
                    fprintf('Processing PDG maps.\n');
                    pdg_cellInspection = obj.visualCellInspection_preprocess(pdg_maps, box_size, cellMasks);

                    obj.Vis_cellInspection_maps = structPacker(obj.Vis_cellInspection_maps,...
                                                                        natmov_cellInspection, 'NatMov', pdg_cellInspection, 'PDG');
                end
                if savedata_flag == 1
                    cellJudging_function(obj.Vis_cellInspection_maps, startcell, cellinfo);
                else
                    cellJudging_function(obj.Vis_cellInspection_maps);
                end
            end
            if strcmp(use_criteria, 'Visual_Inspection_upload')
                %called after step A, to import cellInfo results to roiINFO property
                cellInfo_data = importdata('cellInfo.mat');
                obj.RoiINFO.quality = cellInfo_data.quality;
                obj.RoiINFO.presence = cellInfo_data.presence;
            end
            
            obj.RoiINFO.all = ones(1, obj.num_cells);
            %add more methods here for determining which cells to use in analysis            %Various roi usage criteria
        end
        
        function saveData(obj, varargin)
        % input args are names of the data structures to be saved, or 'all' for all of them
            if any(strcmp(varargin, 'RespData')) || (any(strcmp(varargin, 'all')) && ~isempty(fieldnames(obj.RespData)))
                fprintf('Saving response data...\n');
                RespData = obj.RespData;
                try
                    save('RespData', 'RespData', '-append');
                catch
                    save('RespData', 'RespData');
                end
            end
            
            if any(strcmp(varargin, 'OrientationData')) || (any(strcmp(varargin, 'all')) && ~isempty(fieldnames(obj.OrientationData)))
                fprintf('Saving ori tuning data...\n');
                OrientationData = obj.OrientationData;
                try
                    save('OrientationData', 'OrientationData', '-append');
                catch
                    save('OrientationData', 'OrientationData');
                end
            end

            if any(strcmp(varargin, 'StabilityData')) || (any(strcmp(varargin, 'all')) && ~isempty(fieldnames(obj.StabilityData)))
                fprintf('Saving stability data...\n');
                StabilityData = obj.StabilityData;
                try
                    save('StabilityData', 'StabilityData', '-append');
                catch
                    save('StabilityData', 'StabilityData');
                end
            end

            if any(strcmp(varargin, 'RFdata')) || (any(strcmp(varargin, 'all')) && ~isempty(fieldnames(obj.RFdata)))
                fprintf('Saving receptive field data...\n');
                RFdata = obj.RFdata;
                try
                    save('RFdata', 'RFdata', '-append');
                catch
                    save('RFdata', 'RFdata');
                end
            end

            if any(strcmp(varargin, 'RoiINFO')) || (any(strcmp(varargin, 'all')) && ~isempty(fieldnames(obj.RoiINFO)))
                fprintf('Saving responsiveness info...\n');
                RoiINFO = obj.RoiINFO;
                try
                    save('RoiINFO', 'RoiINFO', '-append');
                catch
                    save('RoiINFO', 'RoiINFO');
                end
            end 
        end
        
        %Getters
        function out = getOrientationData(obj)
            out = obj.OrientationData;
        end
        
        function out = getRespData(obj)
            out = obj.RespData;
        end
        
        function out = getStabilityData(obj)
            out = obj.StabilityData;
        end
        
        function out = getRFdata(obj)
            out = obj.RFdata;
        end
        
        function RoiINFO = getRoiINFO(obj)
            RoiINFO = obj.RoiINFO;
        end

        function setSubsampleFlag(obj, input)
            obj.subsample_flag = input;
        end

        function [reliable, top_percentiles] = trialShuffleTest(obj, Resp, CCs, use_full)
            %trial shuffling to determine responsiveness
            % Input args: 
            %       Resp - input response matrix from RespData 
            %       CCs - real reliability information from StabilityData
            %       use_full - 'Yes' to use full response matrix (all sessions), 'No' to just use D0 (default)
            shuffle_iterations = 1000;
            subsample_iterations = 10;
            percentile = obj.PERCENTILE;
            top_percentiles = zeros(1, obj.num_cells);
            if nargin < 4
                use_full = 'No';            % whether or not to use full respmat or just D0 (D0 for all but stim gap experiments)
            end

            if strcmp(use_full, 'Yes')
                Resp_trialcat = obj.catTrials(Resp);        % full
            else
                Resp_trialcat = squeeze(Resp(:, :, :, 1));       % D0
            end
            if ~strcmp(obj.subsample_flag, 'Yes')
                reliable = zeros(1, obj.num_cells);
                num_reps = size(Resp_trialcat, 1);
                for ii = 1:obj.num_cells
                    ii
                    
                    curr_CCs = zeros(1, shuffle_iterations);
                    actual_CC = obj.btwTrialCC(Resp_trialcat(:, :, ii));  
        
                    for kk = 1:shuffle_iterations
                        shuffled_Resp = zeros(num_reps, size(Resp_trialcat, 2));
                        for rr = 1:num_reps
                            shuffled_Resp(rr, :) = circshift(squeeze(Resp_trialcat(rr, :, ii)), randi(size(Resp_trialcat, 2)), 2);
                        end
                        curr_CCs(kk) = obj.btwTrialCC(shuffled_Resp);
                    end
                    thresh = prctile(curr_CCs, percentile);
                    reliable(ii) = thresh < actual_CC;
                    top_percentiles(ii) = thresh;
                end
            else
                if strcmp(use_full, 'Yes')
                    obj.subsample = obj.findTrialcatSubsample;          % if using full respmat, need to scale subsampling based on number of trials
                else
                    obj.subsample = 8;                                  % if only using D0, subsample 8 to match PDG
                end
                reliable = zeros(subsample_iterations, obj.num_cells);
                %generate list of trial subsamples
                for bb = 1:subsample_iterations
                    trialshuffle = randperm(size(Resp_trialcat, 1));
                    trialselect(bb, :) = trialshuffle(1:obj.subsample);
                end
                %generate actual CC
                actual_CC_temp = zeros(subsample_iterations, obj.num_cells); 
                for bb = 1:subsample_iterations
                    RespMat_subsampled = Resp_trialcat(trialselect(bb, :), :, :);
                    for ii = 1:obj.num_cells
                        actual_CC_temp(bb, ii) = obj.btwTrialCC(RespMat_subsampled(:, :, ii));
                    end
                end
                actual_CC = mean(actual_CC_temp, 1);

                %generate shuffle-trial CC
                for bb = 1:subsample_iterations
                    bb
                    RespMat_subsampled = Resp_trialcat(trialselect(bb, :), :, :);
                    num_reps = size(RespMat_subsampled, 1);
                    for ii = 1:obj.num_cells
                        % ii
                        curr_CCs = zeros(1, shuffle_iterations);
                        for kk = 1:shuffle_iterations
                            shuffled_Resp = zeros(num_reps, size(RespMat_subsampled, 2));
                            for rr = 1:num_reps
                                shuffled_Resp(rr, :) = circshift(squeeze(RespMat_subsampled(rr, :, ii)), randi(size(RespMat_subsampled, 2)), 2);
                            end
                            curr_CCs(kk) = obj.btwTrialCC(shuffled_Resp);
                        end
                        threshold(bb, ii) = prctile(curr_CCs, percentile);
                    end
                end
                avg_threshold = mean(threshold, 1);
                reliable = actual_CC > avg_threshold;  
                top_percentiles = avg_threshold;         
            end
        end

        function out = catTrials(obj, Resp)
            out = [];
            for kk = 1:size(Resp, 4)
                out = cat(1, out, Resp(:, :, :, kk));          
            end
        end
    end

    
    methods (Access = private)
        
        function cell_maps = visualCellInspection_preprocess(obj, maps, box_size, cellMasks)
        % preprocessor for preparing cell images for visual inspection software
            if nargin < 2
                fprintf('Select map file.\n')
                filename = uigetfile('.mat', 'MultiSelect', 'on');
                maps = importdata(filename);
                if nargin < 3
                    box_size = 40;         %box_size x box_size pixels
                    if nargin < 4
                        fprintf('Select reference data file.\n');
                        reference_file = uigetfile('.mat');
                        reference_data = importdata(reference_file);
                        cellMasks = reference_data.cellMasks;
                    end
                end
            end
            avg_projections = maps.avg_projections;
            activity_maps = maps.activity_maps;
            obj.num_sessions = size(avg_projections, 3);
            yPixels = size(avg_projections, 1);
            xPixels = size(avg_projections, 2);
            
            %set default width for cell crop boxes
            padding = 100;
            
            % normalize, store in matrix
            for ii = 1:obj.num_sessions
                curr_avgproj = avg_projections(:, :, ii);
                curr_actmap = activity_maps(:, :, ii);
                if ii == 1
                    avgproj_mat = zeros(yPixels+(padding*2), xPixels+(padding*2), obj.num_sessions);        %initialize padded matrix
                    actmap_mat = zeros(yPixels+(padding*2), xPixels+(padding*2), obj.num_sessions);        %initialize padded matrix
                end
                curr_avgproj = obj.normalize_map(curr_avgproj);     %normalize
                curr_actmap = obj.normalize_map(curr_actmap);
                avgproj_mat(padding+1:padding+yPixels, padding+1:padding+xPixels, ii) = curr_avgproj;
                actmap_mat(padding+1:padding+yPixels, padding+1:padding+xPixels, ii) = curr_actmap;
            end

            %% extract cell array of cropped individual cell maps across all sessions
            obj.num_cells = length(cellMasks);
            cell_avgproj_mat = zeros(box_size, box_size, obj.num_sessions, obj.num_cells);        
            cell_actmap_mat = zeros(box_size, box_size, obj.num_sessions, obj.num_cells);        

            %find centers
            centers = zeros(obj.num_cells, 2);
            for jj = 1:obj.num_cells
                currCenter = obj.findCentroid(jj, cellMasks);
                centers(jj, :) = currCenter+padding;           %adjust coords based on padding
            end
            
            %extract cropped map of each cell across sessions
            fprintf('Extracting cropped cells...\n');
            for kk = 1:obj.num_cells
                y_lower = centers(kk, 2) - box_size/2;          
                y_upper = centers(kk, 2) + box_size/2;
                x_lower = centers(kk, 1) - box_size/2;
                x_upper = centers(kk, 1) + box_size/2;
                for yy = 1:obj.num_sessions
                    temp = obj.normalize_map(avgproj_mat(y_lower+1:y_upper, x_lower+1:x_upper, yy));
                    scaleto = min(temp(temp ~= 0));
                    temp(temp == 0 ) = scaleto;     %scaling padding to match mininum value of cropped box
                    cell_avgproj_mat(:, :, yy, kk) = temp;           
                    cell_actmap_mat(:, :, yy, kk) = obj.normalize_map(actmap_mat(y_lower+1:y_upper, x_lower+1:x_upper, yy));
                end
            end
            
            cell_maps = structPacker([], avgproj_mat, 'avgproj_mat', actmap_mat, 'actmap_mat',...
                                    cell_avgproj_mat, 'cell_avgproj_mat', cell_actmap_mat, 'cell_actmap_mat', centers, 'centers',...
                                    xPixels, 'xPixels', yPixels, 'yPixels', padding, 'padding');
        end

        function map_out = normalize_map(obj, map_in)
            tmp = map_in - min(min(map_in));
            map_out = tmp/max(max(tmp));
        end
        
        function center = findCentroid(obj, cellID, cellMasks)
            roiCoords = cellMasks{cellID};
            center = [round(mean(roiCoords(:, 1))), round(mean(roiCoords(:, 2)))];
        end

        function CC = btwTrialCC(obj, Resp)
            
            num_reps = size(Resp, 1);
            CC_data_iterated = zeros(1, num_reps);
            for rep = 1:num_reps
                CC_data_iterated(rep) = corr(Resp(rep,:)', nanmean(Resp(1:end~=rep,:),1)');
            end
            CC = nanmean(CC_data_iterated);
            
            %basic pearson CC calculation using averaging of 2 sets of random trials (iterated a number of times)
            % num_reps = size(Resp, 1);
            % iterations = 10;
            % CC_data_iterated = zeros(1, iterations);
            % for qq = 1:iterations
            %     trialshuffle = randperm(num_reps);
            %     trialselect1 = trialshuffle(1:floor(num_reps/2));
            %     trialselect2 = trialshuffle(floor(num_reps/2)+1:end);
            %     randavg1 = mean(Resp(trialselect1, :), 1);
            %     randavg2 = mean(Resp(trialselect2, :), 1);
            %     CC_data_iterated(qq) = corr(randavg1', randavg2');
            % end
            % CC = nanmean(CC_data_iterated);
        end

        function out = findTrialcatSubsample(obj)
            out = size(obj.RespData.PDG.RespMat_Full, 1)*obj.num_sessions;
        end


    end
end
                    
                
                