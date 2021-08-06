classdef PDG_Processor < ResponseProcessor & StabilityProcessor
    
    properties (Access = protected)
        OrientationData                     % data structure containing orientation tuning data
        ori_flag = 'Yes';                   % whether or not to run orientation tuning analysis to include in data struct
    end
    
    methods
        function obj = PDG_Processor()
            obj.OrientationData = struct();
        end
        
        function run(obj)
            obj.getFiles('PDG');
            obj.Stimdata = obj.getStimData('PDG', obj.ref_session);
            
            %Sort into Response matrix and calculate btw trial CCs for each session
            RespMat_Full = obj.sortData;            
            obj.RespData = structPacker(obj.RespData, RespMat_Full, 'RespMat_Full');        
         
            %Gather orientation tuning statistics
            if strcmp(obj.ori_flag, 'Yes')
                obj.extractOriTuning();
            end
               
            %Instability Index calculations
            CCs = obj.computeCCs(RespMat_Full);
            [RDI, CC_ws, CC_bs, RDI_control] = obj.computeStability(RespMat_Full);
            obj.StabilityData = structPacker(obj.StabilityData, CCs, 'CCs', RDI, 'RDI', CC_ws, 'CC_ws', CC_bs, 'CC_bs', RDI_control, 'RDI_control');
        end
        
        function extractOriTuning(obj)
            isTuned = zeros(obj.num_sessions, obj.num_cells);
            osiMat = zeros(obj.num_sessions, obj.num_cells);
            zScoreMat = zeros(obj.num_sessions, obj.num_cells);
            oriPref = zeros(obj.num_sessions, obj.num_cells);
            oriResp = zeros(obj.num_sessions, obj.num_cells, obj.Stimdata.setsperRepeat);
            
            for kk = 1:obj.num_sessions
                fprintf('Processing ori tuning: %d out of %d\n', kk, obj.num_sessions);
                if iscell(obj.filelist)
                    filename = obj.filelist{kk};
                else
                    filename = obj.filelist;
                end
                OrientationTuning_stab(filename, 'None', 'Yes');
                currdata = importdata(filename);
                
                isTuned(kk, :) = currdata.isTuned;          %save relevant statistics
                osiMat(kk, :) = currdata.osiVec;
                zScoreMat(kk, :) = currdata.zScoreVec;
                oriResp(kk, :, :) = currdata.oriResp;
                oriPref(kk, :) = currdata.oriPref;
            end
            
            tuned_thresh = obj.num_sessions-2;       %uses cells tuned for tuned_thresh number of sessions out of total
            osiAvg = mean(osiMat(:, sum(isTuned, 1)>=tuned_thresh), 2);
            zScoreAvg = mean(zScoreMat(:, sum(isTuned, 1)>=tuned_thresh), 2);
            osiSEM = std(osiMat(:, sum(isTuned,1)>=tuned_thresh), [], 2)/sqrt(sum(sum(isTuned,1)>=tuned_thresh));
            zScoreSEM = std(zScoreMat(:, sum(isTuned,1)>=tuned_thresh), [], 2)/sqrt(sum(sum(isTuned,1)>=tuned_thresh));
            
            obj.OrientationData = structPacker(obj.OrientationData, isTuned, 'isTuned', osiMat, 'osiMat', zScoreMat, 'zScoreMat',...
                oriResp, 'oriResp', oriPref, 'oriPref', osiAvg, 'osiAvg', zScoreAvg, 'zScoreAvg', osiSEM, 'osiSEM', zScoreSEM, 'zScoreSEM');  
        end
        
        %Getters
        function OrientationData = getOrientationData(obj)
            OrientationData = obj.OrientationData;
        end

        %Setters
        function setOriflag(obj, flag)
            obj.ori_flag = flag;
        end
    end
end