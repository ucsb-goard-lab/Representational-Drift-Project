classdef ResponseProcessor < General_Processor
    
    properties (Access = protected)
        RespData                        % neural response matrix data structure
        dataType                        % type of data in use ('DFF' vs 'spikes')
    end
    
    methods
        
        function obj = ResponseProcessor()
            obj.RespData = struct();
        end

        function RespMat = sortData(obj)
            num_stimfiles = length(obj.Stimdata);
            if num_stimfiles > 1
                RespMat = zeros(obj.Stimdata(1).totalReps, obj.Stimdata(1).repeatTime*obj.Stimdata(1).framerate, obj.num_cells, obj.num_sessions);  %reps, frames, cells, sessions          
            else
                RespMat = zeros(obj.Stimdata.totalReps, obj.Stimdata.repeatTime*obj.Stimdata.framerate, obj.num_cells, obj.num_sessions);  %reps, frames, cells, sessions          
            end
            
            for kk = 1:obj.num_sessions
                fprintf('Importing session %d out of %d\n', kk, obj.num_sessions);
                if iscell(obj.filelist)
                    filename = obj.filelist{kk};
                else
                    filename = obj.filelist;
                end
                currdata = importdata(filename);
                
                if num_stimfiles > 1
                    curr_stimdata = obj.Stimdata(kk);
                else
                    curr_stimdata = obj.Stimdata;
                end
                fprintf('Sorting data...\n')
                RespMat(:, :, :, kk) = obj.responseSorter(currdata, curr_stimdata);
            end
            clearvars currdata
        end
        
        function RespMat = responseSorter(obj, data, stimdata)
            %general repeat sorter

            repeatFrames = round(stimdata.repeatTime*stimdata.framerate);
            totalReps = stimdata.totalReps;

            RespMat = zeros(totalReps, repeatFrames, obj.num_cells);

            for rep = 1:totalReps
                currFrame = (rep-1)*repeatFrames;
                switch obj.dataType
                    case 'DFF'
                        RespMat(rep, :, :) = data.DFF(:, currFrame+1:currFrame+repeatFrames)';
                    case 'spikes'
                        RespMat(rep, :, :) = data.spikes(:, currFrame+1:currFrame+repeatFrames)';
                end
            end

            try                                                 %descramble data if necessary (if trialOrder is provided)
                [~, sorting_vector] = sort(stimdata.trialOrder);
                RespMat = RespMat(sorting_vector, :, :);
            catch
            end
        end 
        
        function onTime_RespMat = extractOnTime(obj, Full_RespMat)
            for kk = 1:obj.num_sessions
                if length(obj.Stimdata) > 1
                    currstimdata = obj.Stimdata(kk);
                else
                    currstimdata = obj.Stimdata;
                end
                onTime_RespMat(:, :, :, kk) = Full_RespMat(:, (currstimdata.offTime*currstimdata.framerate)+1:...
                                                        (currstimdata.offTime*currstimdata.framerate)+...
                                                        (currstimdata.onTime*currstimdata.framerate), :, kk);
            end
        end
        
        %Getter
        function RespData = getRespData(obj)
            RespData = obj.RespData;
        end
        
        %Setter
        function setRespData(obj, input)
            obj.RespData = input;
        end

        function setDataType(obj, input)
            obj.dataType = input;
        end
    end
end