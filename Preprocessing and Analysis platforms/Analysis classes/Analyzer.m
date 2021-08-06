classdef Analyzer < ChronicImaging_Mouse

    properties
        max_sessions        %highest number of sessions of all added data
        num_cells_bymouse      %number of cells in each concatenated mouse (in the order in which they were concatenated)  
        num_sessions_bymouse    % for number of sessions
    end

    properties (Access = protected)
        quality_threshold   %subjective cell quality threshold for cell inclusion
        use_gates            %cell inclusion criteria for last called method
        use_cells            %last usecells vector called
        use_sessions_full         %last use_sessions matrix called
        last_amt_added      %last amount of cells added by concatenator method
        recording_weeks     % [1 x sessions] indicates if a mouse has a recording on a given week
    end

    properties (Access = private)
    end

    methods
        function obj = Analyzer()     
            %empty
        end

        %combining data from multiple windows
        function addData(obj, weeks)           %add ability for input args (as in constructor)
            if isempty(obj.max_sessions)
                obj.max_sessions = obj.num_sessions;
            end
            if isempty(obj.num_cells_bymouse)
                obj.num_cells_bymouse = obj.num_cells;
            end
            if isempty(obj.last_amt_added)
                obj.last_amt_added = obj.num_cells;
            end
            if isempty(obj.num_sessions_bymouse)
                obj.num_sessions_bymouse = obj.num_sessions;
            end

            if nargin < 2
                weeks = [];
            end
 
            proplist =  properties(obj);
            disp('Select all new data files')
            dataNames = uigetfile('.mat', 'MultiSelect', 'on');
            if ~iscell(dataNames)
                [~, loc] = ismember(erase(dataNames, '.mat'), proplist);
                [obj.(proplist{loc}), num_cells_added, num_sessions_added] = obj.concatenateData(dataNames, weeks);         %outputs concatenated fields of each input and the number of cells added on to previous value
            else
                for zz = 1:length(dataNames)
                    [~, loc] = ismember(erase(dataNames{zz}, '.mat'), proplist);
                    [obj.(proplist{loc}), num_cells_added, num_sessions_added] = obj.concatenateData(dataNames{zz}, weeks);
                end
            end
            if num_cells_added ~= obj.last_amt_added                     %check if this window has already been processed by comparing cell count, if not then update values
                obj.num_cells = obj.num_cells + num_cells_added;
                obj.num_cells_bymouse = [obj.num_cells_bymouse num_cells_added];           %store each window's num_cells value
                obj.num_sessions_bymouse = [obj.num_sessions_bymouse num_sessions_added]
            end
            obj.last_amt_added = num_cells_added;
        end
        
        %gating
        function cellSelection(obj)
            obj.setGates(obj.selectGates);      %select logical gates (criteria)
            
            [cells, sessions] = obj.combineGates(obj.getQualityThreshold);     %combine criteria and output logicals 
            
            obj.setUse_cells(cells);
            obj.setUse_sessions(sessions);
        end   
        
        %Getters
        function thresh = getQualityThreshold(obj)
            thresh = obj.quality_threshold;
        end
        
        function use_cells = getUse_cells(obj)
            use_cells = obj.use_cells;
        end
        
        function use_sessions = getUse_sessions(obj)
            use_sessions = obj.use_sessions_full;
        end

        function weeks = getWeeks(obj)
            weeks = obj.recording_weeks;
        end
        
        function cells = getCellList(obj)
            cells = obj.num_cells_bymouse;
        end
        
        %Setters
        function setQualityThreshold(obj, thresh)
            obj.quality_threshold = thresh;
        end
        
        function setGates(obj, gates)
            obj.use_gates = gates;
        end
        
        function setUse_cells(obj, cells)
            obj.use_cells = cells;
        end
        
        function setUse_sessions(obj, sessions)
            obj.use_sessions_full = sessions;
        end

        function setWeeks(obj, weeks)
            obj.recording_weeks = weeks;
        end

        function corrected = alignWeeks(obj, indata, ses_dim, cell_dim, datatype, cf)
            % for mice that have periods of non-consecutive sessions
            % indata: input data matrix
            % ses_dim: session dimension of input data (can prob automate this, but manual input for now)
            % cell_dim: cells dimension of input data
            % datatype: which kind of data struct the input data comes from, e.g. "RespData"
            % cf: if data starts on session 2, input cf = -1, else disregard
            if ~exist('cf', 'var')
                cf = 0;
            end
            try
                num_cells = size(indata, cell_dim);         % number of cells 
            catch
                num_cells = 0;
            end

            shift_idx = find(diff(obj.recording_weeks) == -1)+ 1 + cf;       
            if ~isempty(shift_idx)
                switch datatype
                    case 'StabilityData'
                        insert_vec = NaN(1, num_cells);             % going to insert a vector of NaNs in the shift index location, and move everything to the right
                        if ses_dim == 1
                            corrected = cat(1,  indata(1:shift_idx-1, :), insert_vec, indata(shift_idx:end, :));
                        elseif ses_dim == 2
                            corrected = cat(2, indata(:, 1:shift_idx-1), insert_vec, indata(:, shift_idx:end));
                        end
                    case 'OrientationData'
                        if num_cells ~= 0
                            if length(size(indata)) == 3
                                insert_vec = NaN(1, num_cells, size(indata, 3));
                                corrected = cat(1,  indata(1:shift_idx-1, :, :), insert_vec, indata(shift_idx:end, :, :));
                            else
                                insert_vec = NaN(1, num_cells);
                                corrected = cat(1,  indata(1:shift_idx-1, :), insert_vec, indata(shift_idx:end, :));
                            end
                        else
                            insert_vec = NaN;
                            corrected = cat(1, indata(1:shift_idx-1), insert_vec, indata(shift_idx:end));
                        end
                    case 'RespData'
                        insert_vec = NaN(size(indata, 1), size(indata, 2), num_cells);
                        corrected = cat(4, indata(:, :, :, 1:shift_idx-1), insert_vec, indata(:, :, :, shift_idx:end));
                    case 'RFdata'
                        if num_cells == 0
                            insert_vec = NaN;
                            corrected = cat(2,  indata(1:shift_idx-1), insert_vec, indata(shift_idx:end));
                        elseif ses_dim == 5
                            insert_vec = NaN(num_cells, size(indata, 2), size(indata, 3), size(indata, 4));
                            corrected = cat(5, indata(:, :, :, :, 1:shift_idx-1), insert_vec, indata(:, :, :, :, shift_idx:end));
                        else
                            if ses_dim == 1
                                corrected = cat(1,  indata(1:shift_idx-1, :), insert_vec, indata(shift_idx:end, :));
                            elseif ses_dim == 2
                                corrected = cat(2, indata(:, 1:shift_idx-1), insert_vec, indata(:, shift_idx:end));
                            end
                        end
                    case 'RoiINFO'
                        insert_vec = NaN(1, num_cells)';             % going to insert a vector of NaNs in the shift index location, and move everything to the right
                        corrected = cat(2, indata(:, 1:shift_idx-1), insert_vec, indata(:, shift_idx:end));
                    case 'pupilinfo'
                        insert_vec = struct();          %empty struct
                        insert_vec.raw = [];
                        insert_vec.resampled = [];
                        insert_vec.sorted = [];
                        corrected = cat(2, indata(1:shift_idx-1), insert_vec, indata(shift_idx:end));
                end
            else
                corrected = indata;
            end
        end  
    end
        
    methods (Access = private)
        
        function [out, num_cells_added, num_sessions_added] = concatenateData(obj, new_filename, weeks)       %move this to private when finished
            dataType = erase(new_filename, '.mat');
            currData = importdata(new_filename);
            PDG_flag = any(strcmp(fieldnames(currData), 'PDG'));
            NatMov_flag = any(strcmp(fieldnames(currData), 'NatMov'));
            PDGcont_flag = any(strcmp(fieldnames(currData), 'PDGcont'));
            NatMovdisc_flag = any(strcmp(fieldnames(currData), 'NatMovdisc'));
            MultiMov_flag = any(strcmp(fieldnames(currData), 'MultiMov'));
            switch dataType
                case 'RespData'
                    fprintf('Adding response data...\n');
                    num_cells_added = size(currData.PDG.RespMat_Full, 3);
                    num_sessions_added = size(currData.PDG.RespMat_Full, 4);
                   
                    if ~isempty(weeks)
                        obj.setWeeks(weeks);
                        if PDG_flag
                            currData.PDG.RespMat_Full = obj.alignWeeks(currData.PDG.RespMat_Full, 4, 3, 'RespData');
                        end
                        if NatMov_flag
                            currData.NatMov.RespMat_Full = obj.alignWeeks(currData.NatMov.RespMat_Full, 4, 3, 'RespData');
                            currData.NatMov.RespMat_onTime = obj.alignWeeks(currData.NatMov.RespMat_onTime, 4, 3, 'RespData');
                        end
                        if PDGcont_flag
                            currData.PDGcont.RespMat_Full = obj.alignWeeks(currData.PDGcont.RespMat_Full, 4, 3, 'RespData');
                            currData.PDGcont.RespMat_onTime = obj.alignWeeks(currData.PDGcont.RespMat_onTime, 4, 3, 'RespData');
                        end
                        if NatMovdisc_flag
                            currData.NatMovdisc.RespMat_Full = obj.alignWeeks(currData.NatMovdisc.RespMat_Full, 4, 3, 'RespData');
                        end
                    end



                    if PDG_flag
                        curr_sessions = size(currData.PDG.RespMat_Full, 4);
                    elseif NatMov_flag
                        curr_sessions = size(currData.NatMov.RespMat_Full, 4);
                    elseif MultiMov_flag
                        curr_sessions = size(currData.MultiMov.RespMat_Full{1}, 4);
                    end

                    if curr_sessions > obj.max_sessions
                        obj.max_sessions = curr_sessions;
                    end

                    if curr_sessions < obj.max_sessions

                        if PDG_flag
                            currData.PDG.RespMat_Full(:, :, :, curr_sessions+1:obj.max_sessions) = NaN;
                        end

                        if NatMov_flag
                            currData.NatMov.RespMat_Full(:, :, :, curr_sessions+1:obj.max_sessions) = NaN;
                            currData.NatMov.RespMat_onTime(:, :, :, curr_sessions+1:obj.max_sessions) = NaN; 
                        end

                        if PDGcont_flag
                            currData.PDGcont.RespMat_Full(:, :, :, curr_sessions+1:obj.max_sessions) = NaN;
                            currData.PDGcont.RespMat_onTime(:, :, :, curr_sessions+1:obj.max_sessions) = NaN; 
                        end

                        if NatMovdisc_flag
                            currData.NatMovdisc.RespMat_Full(:, :, :, curr_sessions+1:obj.max_sessions) = NaN;
                        end

                        if MultiMov_flag
                            for mm = 1:length(currData.MultiMov.RespMat_Full)
                                currData.MultiMov.RespMat_Full{mm}(:, :, :, curr_sessions+1:obj.max_sessions) = NaN;
                                currData.MultiMov.RespMat_onTime{mm}(:, :, :, curr_sessions+1:obj.max_sessions) = NaN; 
                            end
                        end
                                    %%% add more stuff here for RF
                    end
                       
                    if PDG_flag
                        out.PDG.RespMat_Full = cat(3, obj.RespData.PDG.RespMat_Full, currData.PDG.RespMat_Full);
                    end
                    if NatMov_flag
                        if size(currData.NatMov.RespMat_Full, 1) == size(obj.RespData.NatMov.RespMat_Full, 1)
                            out.NatMov.RespMat_Full = cat(3, obj.RespData.NatMov.RespMat_Full, currData.NatMov.RespMat_Full);
                            out.NatMov.RespMat_onTime = cat(3, obj.RespData.NatMov.RespMat_onTime, currData.NatMov.RespMat_onTime);
                        else            % this accounts for adding data with fewer repeats than normal
                            dimdiff = size(obj.RespData.NatMov.RespMat_Full, 1) - size(currData.NatMov.RespMat_Full, 1);
                            currData.NatMov.RespMat_Full = cat(1, currData.NatMov.RespMat_Full, nan(dimdiff, size(currData.NatMov.RespMat_Full, 2), size(currData.NatMov.RespMat_Full, 3), size(currData.NatMov.RespMat_Full, 4)));
                            currData.NatMov.RespMat_onTime = cat(1, currData.NatMov.RespMat_onTime, nan(dimdiff, size(currData.NatMov.RespMat_onTime, 2), size(currData.NatMov.RespMat_onTime, 3), size(currData.NatMov.RespMat_onTime, 4)));
                            out.NatMov.RespMat_Full = cat(3, obj.RespData.NatMov.RespMat_Full, currData.NatMov.RespMat_Full);
                            out.NatMov.RespMat_onTime = cat(3, obj.RespData.NatMov.RespMat_onTime, currData.NatMov.RespMat_onTime);
                        end
                    end
                    if PDGcont_flag
                        out.PDGcont.RespMat_Full = cat(3, obj.RespData.PDGcont.RespMat_Full, currData.PDGcont.RespMat_Full);
                        out.PDGcont.RespMat_onTime = cat(3, obj.RespData.PDGcont.RespMat_onTime, currData.PDGcont.RespMat_onTime);
                    end
                    if NatMovdisc_flag
                        out.NatMovdisc.RespMat_Full = cat(3, obj.RespData.NatMovdisc.RespMat_Full, currData.NatMovdisc.RespMat_Full);
                    end
                    if MultiMov_flag
                        for mm = 1:length(currData.MultiMov.RespMat_Full)
                            out.MultiMov.RespMat_Full{mm} = cat(3, obj.RespData.MultiMov.RespMat_Full{mm}, currData.MultiMov.RespMat_Full{mm});
                            out.MultiMov.RespMat_onTime{mm} = cat(3, obj.RespData.MultiMov.RespMat_onTime{mm}, currData.MultiMov.RespMat_onTime{mm});
                        end
                    end
                    %%% add more stuff here for RF
                case 'StabilityData'
                    fprintf('Adding stability data...\n');
                    num_cells_added = size(currData.PDG.CCs, 2);
                    num_sessions_added = size(currData.PDG.CCs, 1);

                    if ~isempty(weeks)
                        obj.setWeeks(weeks);
                        if PDG_flag
                            currData.PDG.RDI = obj.alignWeeks(currData.PDG.RDI, 1, 2, 'StabilityData', -1);           % correction factor of -1 due to first index being week 2
                            currData.PDG.CCs = obj.alignWeeks(currData.PDG.CCs, 1, 2, 'StabilityData');
                            currData.PDG.CC_ws = obj.alignWeeks(currData.PDG.CC_ws, 1, 2, 'StabilityData', -1);
                            currData.PDG.CC_bs = obj.alignWeeks(currData.PDG.CC_bs, 1, 2, 'StabilityData', -1);
                        end
                        if NatMov_flag
                            currData.NatMov.RDI = obj.alignWeeks(currData.NatMov.RDI, 1, 2, 'StabilityData', -1);           % correction factor of -1 due to first index being week 2
                            currData.NatMov.CCs = obj.alignWeeks(currData.NatMov.CCs, 1, 2, 'StabilityData');
                            currData.NatMov.CC_ws = obj.alignWeeks(currData.NatMov.CC_ws, 1, 2, 'StabilityData', -1);
                            currData.NatMov.CC_bs = obj.alignWeeks(currData.NatMov.CC_bs, 1, 2, 'StabilityData', -1);                        
                        end
                        if PDGcont_flag
                            currData.PDGcont_flag.RDI = obj.alignWeeks(currData.PDGcont_flag.RDI, 1, 2, 'StabilityData', -1);           % correction factor of -1 due to first index being week 2
                            currData.PDGcont_flag.CCs = obj.alignWeeks(currData.PDGcont_flag.CCs, 1, 2, 'StabilityData');
                            currData.PDGcont_flag.CC_ws = obj.alignWeeks(currData.PDGcont_flag.CC_ws, 1, 2, 'StabilityData', -1);
                            currData.PDGcont_flag.CC_bs = obj.alignWeeks(currData.PDGcont_flag.CC_bs, 1, 2, 'StabilityData', -1);     
                        end
                        if NatMovdisc_flag
                            currData.NatMovdisc.RDI = obj.alignWeeks(currData.NatMovdisc.RDI, 1, 2, 'StabilityData', -1);           % correction factor of -1 due to first index being week 2
                            currData.NatMovdisc.CCs = obj.alignWeeks(currData.NatMovdisc.CCs, 1, 2, 'StabilityData');
                            currData.NatMovdisc.CC_ws = obj.alignWeeks(currData.NatMovdisc.CC_ws, 1, 2, 'StabilityData', -1);
                            currData.NatMovdisc.CC_bs = obj.alignWeeks(currData.NatMovdisc.CC_bs, 1, 2, 'StabilityData', -1); 
                        end
                    end 

                    if PDG_flag
                        curr_sessions = size(currData.PDG.CCs, 1);
                    elseif NatMov_flag
                        curr_sessions = size(currData.NatMov.CCs, 1);
                    elseif MultiMov_flag
                        curr_sessions = size(currData.MultiMov.CCs{1}, 1);
                    end

                    if curr_sessions > obj.max_sessions
                        obj.max_sessions = curr_sessions;
                    end

                    if curr_sessions < obj.max_sessions 
                        if PDG_flag
                            currData.PDG.RDI(curr_sessions:obj.max_sessions-1, :) = NaN;
                            currData.PDG.CCs(curr_sessions+1:obj.max_sessions, :) = NaN;
                            currData.PDG.CC_ws(curr_sessions:obj.max_sessions-1, :) = NaN;
                            currData.PDG.CC_bs(curr_sessions:obj.max_sessions-1, :) = NaN;
                        end
                        if NatMov_flag
                            currData.NatMov.RDI(curr_sessions:obj.max_sessions-1, :) = NaN;
                            currData.NatMov.CCs(curr_sessions+1:obj.max_sessions, :) = NaN; 
                            currData.NatMov.CC_ws(curr_sessions:obj.max_sessions-1, :) = NaN; 
                            currData.NatMov.CC_bs(curr_sessions:obj.max_sessions-1, :) = NaN;
                        end
                        if PDGcont_flag
                            currData.PDGcont.RDI(curr_sessions:obj.max_sessions-1, :) = NaN;
                            currData.PDGcont.CCs(curr_sessions+1:obj.max_sessions, :) = NaN; 
                            currData.PDGcont.CC_ws(curr_sessions:obj.max_sessions-1, :) = NaN; 
                            currData.PDGcont.CC_bs(curr_sessions:obj.max_sessions-1, :) = NaN;
                        end
                        if NatMovdisc_flag
                            currData.NatMovdisc.RDI(curr_sessions:obj.max_sessions-1, :) = NaN;
                            currData.NatMovdisc.CCs(curr_sessions+1:obj.max_sessions, :) = NaN; 
                            currData.NatMovdisc.CC_ws(curr_sessions:obj.max_sessions-1, :) = NaN; 
                            currData.NatMovdisc.CC_bs(curr_sessions:obj.max_sessions-1, :) = NaN;
                        end
                        if MultiMov_flag
                            for mm = 1:length(currData.MultiMov.CCs{1})
                                currData.MultiMov.RDI{mm}(curr_sessions:obj.max_sessions-1, :) = NaN;
                                currData.MultiMov.CCs{mm}(curr_sessions+1:obj.max_sessions, :) = NaN; 
                                currData.MultiMov.CC_ws{mm}(curr_sessions:obj.max_sessions-1, :) = NaN; 
                                currData.MultiMov.CC_bs{mm}(curr_sessions:obj.max_sessions-1, :) = NaN;
                            end
                        end
                    end

                    if PDG_flag
                        out.PDG.RDI = cat(2, obj.StabilityData.PDG.RDI, currData.PDG.RDI);
                        out.PDG.CCs = cat(2, obj.StabilityData.PDG.CCs, currData.PDG.CCs);
                        out.PDG.CC_ws = cat(2, obj.StabilityData.PDG.CC_ws, currData.PDG.CC_ws);
                        out.PDG.CC_bs = cat(2, obj.StabilityData.PDG.CC_bs, currData.PDG.CC_bs);
                        out.PDG.RDI_control = cat(2, obj.StabilityData.PDG.RDI_control, currData.PDG.RDI_control);
                    end
                    if NatMov_flag
                        out.NatMov.RDI = cat(2, obj.StabilityData.NatMov.RDI, currData.NatMov.RDI);
                        out.NatMov.CCs = cat(2, obj.StabilityData.NatMov.CCs, currData.NatMov.CCs);
                        out.NatMov.CC_ws = cat(2, obj.StabilityData.NatMov.CC_ws, currData.NatMov.CC_ws);
                        out.NatMov.CC_bs = cat(2, obj.StabilityData.NatMov.CC_bs, currData.NatMov.CC_bs);
                        out.NatMov.RDI_control = cat(1, obj.StabilityData.NatMov.RDI_control, currData.NatMov.RDI_control);
                    end
                    if PDGcont_flag
                        out.PDGcont.RDI = cat(2, obj.StabilityData.PDGcont.RDI, currData.PDGcont.RDI);
                        out.PDGcont.CCs = cat(2, obj.StabilityData.PDGcont.CCs, currData.PDGcont.CCs);
                        out.PDGcont.CC_ws = cat(2, obj.StabilityData.PDGcont.CC_ws, currData.PDGcont.CC_ws);
                        out.PDGcont.CC_bs = cat(2, obj.StabilityData.PDGcont.CC_bs, currData.PDGcont.CC_bs);
                        out.PDGcont.RDI_control = cat(1, obj.StabilityData.PDGcont.RDI_control, currData.PDGcont.RDI_control);
                    end
                    if NatMovdisc_flag
                        out.NatMovdisc.RDI = cat(2, obj.StabilityData.NatMovdisc.RDI, currData.NatMovdisc.RDI);
                        out.NatMovdisc.CCs = cat(2, obj.StabilityData.NatMovdisc.CCs, currData.NatMovdisc.CCs);
                        out.NatMovdisc.CC_ws = cat(2, obj.StabilityData.NatMovdisc.CC_ws, currData.NatMovdisc.CC_ws);
                        out.NatMovdisc.CC_bs = cat(2, obj.StabilityData.NatMovdisc.CC_bs, currData.NatMovdisc.CC_bs);
                        out.NatMovdisc.RDI_control = cat(2, obj.StabilityData.NatMovdisc.RDI_control, currData.NatMovdisc.RDI_control);
                    end
                    if MultiMov_flag
                        for mm = 1:length(currData.MultiMov.CCs)
                            out.MultiMov.RDI{mm} = cat(2, obj.StabilityData.MultiMov.RDI{mm}, currData.MultiMov.RDI{mm});
                            out.MultiMov.CCs{mm} = cat(2, obj.StabilityData.MultiMov.CCs{mm}, currData.MultiMov.CCs{mm});
                            out.MultiMov.CC_ws{mm} = cat(2, obj.StabilityData.MultiMov.CC_ws{mm}, currData.MultiMov.CC_ws{mm});
                            out.MultiMov.CC_bs{mm} = cat(2, obj.StabilityData.MultiMov.CC_bs{mm}, currData.MultiMov.CC_bs{mm});
                            out.MultiMov.RDI_control{mm} = cat(1, obj.StabilityData.MultiMov.RDI_control{mm}, currData.MultiMov.RDI_control{mm});
                        end
                    end
                case 'OrientationData'
                    fprintf('Adding orientation tuning data...\n');
                    num_cells_added = size(currData.isTuned, 2);
                    num_sessions_added = size(currData.isTuned, 1);
                   
                     if ~isempty(weeks)
                        obj.setWeeks(weeks);
                        currData.isTuned = obj.alignWeeks(currData.isTuned, 1, 2, 'OrientationData');     
                        currData.osiMat = obj.alignWeeks(currData.osiMat, 1, 2, 'OrientationData'); 
                        currData.zScoreMat = obj.alignWeeks(currData.zScoreMat, 1, 2, 'OrientationData');
                        currData.oriResp = obj.alignWeeks(currData.oriResp, 1, 2, 'OrientationData');
                        currData.oriPref = obj.alignWeeks(currData.oriPref, 1, 2, 'OrientationData');
                        currData.osiAvg = obj.alignWeeks(currData.osiAvg, 2, 0, 'OrientationData');
                        currData.zScoreAvg = obj.alignWeeks(currData.zScoreAvg, 2, 0, 'OrientationData');
                        currData.osiSEM = obj.alignWeeks(currData.osiSEM, 2, 0, 'OrientationData');
                        currData.zScoreSEM = obj.alignWeeks(currData.zScoreSEM, 2, 0, 'OrientationData');          
                    end 

                    curr_sessions = size(currData.isTuned, 1);
                    if curr_sessions > obj.max_sessions
                        obj.max_sessions = curr_sessions;
                    end

                    if curr_sessions < obj.max_sessions
                        currData.isTuned(curr_sessions+1:obj.max_sessions, :) = NaN;
                        currData.osiMat(curr_sessions+1:obj.max_sessions, :) = NaN;  
                        currData.zScoreMat(curr_sessions+1:obj.max_sessions, :) = NaN; 
                        currData.oriResp(curr_sessions+1:obj.max_sessions, :, :) = NaN; 
                        currData.oriPref(curr_sessions+1:obj.max_sessions, :) = NaN; 
                        currData.osiAvg(curr_sessions+1:obj.max_sessions) = NaN; 
                        currData.zScoreAvg(curr_sessions+1:obj.max_sessions) = NaN; 
                        currData.osiSEM(curr_sessions+1:obj.max_sessions) = NaN; 
                        currData.zScoreSEM(curr_sessions+1:obj.max_sessions) = NaN; 
                    end

                    out.isTuned = cat(2, obj.OrientationData.isTuned, currData.isTuned);
                    out.osiMat = cat(2, obj.OrientationData.osiMat, currData.osiMat);
                    out.zScoreMat = cat(2, obj.OrientationData.zScoreMat, currData.zScoreMat);
                    out.oriResp = cat(2, obj.OrientationData.oriResp, currData.oriResp);
                    out.oriPref = cat(2, obj.OrientationData.oriPref, currData.oriPref);
                    out.osiAvg = cat(2, obj.OrientationData.osiAvg, currData.osiAvg);
                    out.zScoreAvg = cat(2, obj.OrientationData.zScoreAvg, currData.zScoreAvg);
                    out.osiSEM = cat(2, obj.OrientationData.osiSEM, currData.osiSEM);
                    out.zScoreSEM = cat(2, obj.OrientationData.zScoreSEM, currData.zScoreSEM);
                case 'RFdata'
                    fprintf('Adding receptive field data...\n');
                    num_cells_added = size(currData.RFmapping_results.isSpatiallyTuned, 2);
                    num_sessions_added = size(currData.RFmapping_results.isSpatiallyTuned, 1);
                   
                    if ~isempty(weeks)
                        obj.setWeeks(weeks);
                        currData.RFmapping_results.alt_fit = obj.alignWeeks(currData.RFmapping_results.alt_fit, 1, 2, 'RFdata');     
                        currData.RFmapping_results.azi_fit = obj.alignWeeks(currData.RFmapping_results.azi_fit, 1, 2, 'RFdata'); 
                        currData.RFmapping_results.alt_pref = obj.alignWeeks(currData.RFmapping_results.alt_pref, 1, 2, 'RFdata');
                        currData.RFmapping_results.azi_pref = obj.alignWeeks(currData.RFmapping_results.azi_pref, 1, 2, 'RFdata');
                        currData.RFmapping_results.alt_p = obj.alignWeeks(currData.RFmapping_results.alt_p, 1, 2, 'RFdata');
                        currData.RFmapping_results.azi_p = obj.alignWeeks(currData.RFmapping_results.azi_p, 1, 2, 'RFdata');
                        currData.RFmapping_results.isSpatiallyTuned = obj.alignWeeks(currData.RFmapping_results.isSpatiallyTuned, 1, 2, 'RFdata');
                        currData.RFmapping_results.roi_centroids = obj.alignWeeks(currData.RFmapping_results.roi_centroids, 3, 1, 'RFdata');
                        currData.RFmapping_results.alt_corr = obj.alignWeeks(currData.RFmapping_results.alt_corr, 2, 0, 'RFdata');
                        currData.RFmapping_results.azi_corr = obj.alignWeeks(currData.RFmapping_results.azi_corr, 2, 0, 'RFdata');
                        currData.RFmapping_results.isVisuallyResponsive_alt = obj.alignWeeks(currData.RFmapping_results.isVisuallyResponsive_alt, 2, 1, 'RFdata');
                        currData.RFmapping_results.isVisuallyResponsive_azi = obj.alignWeeks(currData.RFmapping_results.isVisuallyResponsive_azi, 2, 1, 'RFdata');
                        currData.RFmapping_results.RespMat_alt = obj.alignWeeks(currData.RFmapping_results.RespMat_alt, 5, 1, 'RFdata');
                        currData.RFmapping_results.RespMat_azi = obj.alignWeeks(currData.RFmapping_results.RespMat_azi, 5, 1, 'RFdata');
                    end 

                    curr_sessions = size(currData.RFmapping_results.isSpatiallyTuned, 1);
                    if curr_sessions > obj.max_sessions
                        obj.max_sessions = curr_sessions;
                    end

                    if curr_sessions > obj.max_sessions
                        currData.RFmapping_results.alt_fit(curr_sessions+1:obj.max_sessions, :) = NaN;
                        currData.RFmapping_results.azi_fit(curr_sessions+1:obj.max_sessions, :) = NaN;
                        currData.RFmapping_results.alt_pref(curr_sessions+1:obj.max_sessions, :) = NaN;
                        currData.RFmapping_results.azi_pref(curr_sessions+1:obj.max_sessions, :) = NaN;
                        currData.RFmapping_results.alt_p(curr_sessions+1:obj.max_sessions, :) = NaN;
                        currData.RFmapping_results.azi_p(curr_sessions+1:obj.max_sessions, :) = NaN;
                        currData.RFmapping_results.isSpatiallyTuned(curr_sessions+1:obj.max_sessions, :) = NaN;
                        currData.RFmapping_results.roi_centroids(:, :, curr_sessions+1:obj.max_sessions) = NaN;
                        currData.RFmapping_results.alt_corr(curr_sessions+1:obj.max_sessions) = NaN;
                        currData.RFmapping_results.azi_corr(curr_sessions+1:obj.max_sessions) = NaN;
                        currData.RFmapping_results.isVisuallyResponsive_alt(:, curr_sessions+1:obj.max_sessions) = NaN;
                        currData.RFmapping_results.isVisuallyResponsive_azi(:, curr_sessions+1:obj.max_sessions) = NaN;
                        currData.RFmapping_results.RespMat_alt(:, :, :, :, curr_sessions+1:obj.max_sessions) = NaN;
                        currData.RFmapping_results.RespMat_azi(:, :, :, :, curr_sessions+1:obj.max_sessions) = NaN;
                    end

                    out.RFmapping_results.alt_fit = cat(2, obj.RFdata.RFmapping_results.alt_fit, currData.RFmapping_results.alt_fit);
                    out.RFmapping_results.azi_fit = cat(2, obj.RFdata.RFmapping_results.azi_fit, currData.RFmapping_results.azi_fit);
                    out.RFmapping_results.alt_pref = cat(2, obj.RFdata.RFmapping_results.alt_pref, currData.RFmapping_results.alt_pref);
                    out.RFmapping_results.azi_pref = cat(2, obj.RFdata.RFmapping_results.azi_pref, currData.RFmapping_results.azi_pref);
                    out.RFmapping_results.alt_p = cat(2, obj.RFdata.RFmapping_results.alt_p, currData.RFmapping_results.alt_p);
                    out.RFmapping_results.azi_p = cat(2, obj.RFdata.RFmapping_results.azi_p, currData.RFmapping_results.azi_p);
                    out.RFmapping_results.isSpatiallyTuned = cat(2, obj.RFdata.RFmapping_results.isSpatiallyTuned, currData.RFmapping_results.isSpatiallyTuned);
                    out.RFmapping_results.roi_centroids = cat(1, obj.RFdata.RFmapping_results.roi_centroids, currData.RFmapping_results.roi_centroids);
                    out.RFmapping_results.alt_corr = cat(1, obj.RFdata.RFmapping_results.alt_corr, currData.RFmapping_results.alt_corr);
                    out.RFmapping_results.azi_corr = cat(1, obj.RFdata.RFmapping_results.azi_corr, currData.RFmapping_results.azi_corr);
                    out.RFmapping_results.isVisuallyResponsive_alt = cat(1, obj.RFdata.RFmapping_results.isVisuallyResponsive_alt, currData.RFmapping_results.isVisuallyResponsive_alt);
                    out.RFmapping_results.isVisuallyResponsive_azi = cat(1, obj.RFdata.RFmapping_results.isVisuallyResponsive_azi, currData.RFmapping_results.isVisuallyResponsive_azi);
                    out.RFmapping_results.RespMat_alt = cat(1, obj.RFdata.RFmapping_results.RespMat_alt, currData.RFmapping_results.RespMat_alt);
                    out.RFmapping_results.RespMat_azi = cat(1, obj.RFdata.RFmapping_results.RespMat_azi, currData.RFmapping_results.RespMat_azi);
                case 'RoiINFO'
                    fprintf('Adding ROI info...\n');
                    num_cells_added = size(currData.presence, 1);
                    num_sessions_added = size(currData.presence, 2);
                    PDG_flag = any(strncmp(fieldnames(currData), 'PDG_Responsive', 14));
                    NatMov_flag = any(strncmp(fieldnames(currData), 'NatMov_Responsive', 17));
                    MultiMov_flag = any(strncmp(fieldnames(currData), 'MultiMov_Responsive', 19));
                    if ~isempty(weeks)
                        obj.setWeeks(weeks);
                        currData.presence = obj.alignWeeks(currData.presence, 2, 1, 'RoiINFO');     
                    end 

                    curr_sessions = size(currData.presence, 2);
                    if curr_sessions > obj.max_sessions
                        obj.max_sessions = curr_sessions;
                    end

                    if curr_sessions < obj.max_sessions
                        currData.presence(:, curr_sessions+1:obj.max_sessions) = NaN;
                    end

                    if PDG_flag
                        out.PDG_Responsive_shuffle = cat(2, obj.RoiINFO.PDG_Responsive_shuffle, currData.PDG_Responsive_shuffle);
                        out.PDG_Responsive_thresh = cat(2, obj.RoiINFO.PDG_Responsive_thresh, currData.PDG_Responsive_thresh);
                        try
                        out.PDG_Responsive_ontime = cat(2, obj.RoiINFO.PDG_Responsive_ontime, currData.PDG_Responsive_ontime);
                        catch
                        end
                    end
                    if NatMov_flag
                        out.NatMov_Responsive_shuffle = cat(2, obj.RoiINFO.NatMov_Responsive_shuffle, currData.NatMov_Responsive_shuffle);
                        out.NatMov_Responsive_thresh = cat(2, obj.RoiINFO.NatMov_Responsive_thresh, currData.NatMov_Responsive_thresh);
                        try
                        out.NatMov_Responsive_ontime = cat(2, obj.RoiINFO.NatMov_Responsive_ontime, currData.NatMov_Responsive_ontime);
                        catch
                        end
                    end
                    if MultiMov_flag
                        out.MultiMov_Responsive_thresh = cat(2, obj.RoiINFO.MultiMov_Responsive_thresh, currData.MultiMov_Responsive_thresh);
                        out.MultiMov_Responsive_shuffle = cat(2, obj.RoiINFO.MultiMov_Responsive_shuffle, currData.MultiMov_Responsive_shuffle);
                    end

                    out.quality = cat(2, obj.RoiINFO.quality, currData.quality);
                    out.presence = cat(1, obj.RoiINFO.presence, currData.presence);
            end   
        end
        
        function useGates = selectGates(obj)
            possible_criteria = fieldnames(obj.RoiINFO);
            [useGates, ~] = listdlg('ListString', possible_criteria);           
        end
        
        function [useCells, useSessions] = combineGates(obj, quality_threshold)
            presence_flag = 0;
            useCells = ones(1, obj.num_cells);
            if nargin < 2
                if any(strcmp(obj.usegates, 'quality'))
                    quality_threshold = 3;      %default 
                end
            end
            possible_criteria = fieldnames(obj.RoiINFO);
            for ii = 1:length(obj.use_gates)
                currgate_temp = possible_criteria(obj.use_gates(ii));     
                currgate = currgate_temp{1};
                if strcmp(currgate, 'quality')          %if cell quality is included as a criterion, convert quality to logical using threshold
                    quality_logical = obj.RoiINFO.quality >= quality_threshold;
                    useCells = useCells & quality_logical;
                elseif strcmp(currgate, 'presence')
                    presence_flag = 1;
                else
                    useCells = useCells & obj.RoiINFO.(currgate);
                end
            end

            if presence_flag == 1
                useSessions = obj.RoiINFO.presence;
            else
                useSessions = ones(obj.num_cells, obj.num_sessions);      %if presence is not included as a criterion, output ones matrix (use all sessions)
            end
        end

    end
end                
