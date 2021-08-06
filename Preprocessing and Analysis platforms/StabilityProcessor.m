classdef StabilityProcessor < General_Processor
    
    properties (Access = protected)
        CC_method string                    % method for computing reliability
        subsample = 8;                      % default trial subsampling value (matching MOV to PDG)
        subsample_flag;                     % whether or not to subsample trials
        rectify_flag;                       % whether or not to rectify reliability values ('defaults to Yes')
        CC_FUDGE = 0.001;                   % reliability rectification value
    end

    properties
        StabilityData                       % data structure for reliability and stability info
    end
    
    methods
        
        function obj = StabilityProcessor()
            obj.StabilityData = struct();
            obj.rectify_flag = 'Yes';
        end
        
        function CCs = computeCCs(obj, RespMat)

            if ~strcmp(obj.subsample_flag, 'Yes')
                CCs = zeros(obj.num_sessions, obj.num_cells);
                for kk = 1:obj.num_sessions
                    fprintf('Calculating CCs...\n')
                    CCs(kk, :) = obj.Reliability_function(RespMat(:, :, :, kk));
                end
            else
                iterations = 10;
                CCs_temp = zeros(obj.num_sessions, obj.num_cells, iterations);
                for rr = 1:iterations
                    trialshuffle = randperm(size(RespMat, 1));
                    trialselect = trialshuffle(1:obj.subsample);
                    RespMat_subsampled = RespMat(trialselect, :, :, :);
                    for kk = 1:obj.num_sessions
                        fprintf('Calculating CCs...\n')
                        CCs_temp(kk, :, rr) = obj.Reliability_function(RespMat_subsampled(:, :, :, kk));
                    end
                end
                CCs = nanmean(CCs_temp, 3);
            end

            if strcmp(obj.rectify_flag, 'Yes')
                CCs(CCs <= 0) = obj.CC_FUDGE;
            end
        end
        
        function CC_data = Reliability_function(obj, RespMat)
            %CC calculation function
            numReps = size(RespMat, 1);
            CC_data = zeros(1,obj.num_cells);

            switch obj.CC_method
                case 'Average'
                    for cell = 1:obj.num_cells
                        CC_data_curr = zeros(1,numReps);
                        for rep = 1:numReps
                            CC_data_curr(rep) = corr(RespMat(rep,:,cell)',nanmean(RespMat(1:end~=rep,:,cell),1)');
                        end
                        CC_data(cell) = nanmean(CC_data_curr);
                    end
                case 'Evenodd'
                    for cell = 1:obj.num_cells
                        Eventrace = nanmean(RespMat(2:2:end, :, cell), 1);
                        Oddtrace = nanmean(RespMat(1:2:end, :, cell), 1);
                        CC_data(cell) = corr(Eventrace',Oddtrace');
                    end
                case 'Random'
                    if numReps > 8
                        iterations = 100;
                    else
                        iterations = 50;
                    end
                    for cell = 1:obj.num_cells
                        CC_data_iterated = zeros(1, iterations);
                        for jj = 1:iterations
                            trialshuffle = randperm(numReps);
                            trialselect1 = trialshuffle(1:floor(numReps/2));
                            trialselect2 = trialshuffle(floor(numReps/2)+1:end);
                            randavg1 = mean(RespMat(trialselect1, :, cell), 1);
                            randavg2 = mean(RespMat(trialselect2, :, cell), 1);
                            CC_data_iterated(jj) = corr(randavg1', randavg2');
                        end
                        CC_data(cell) = nanmean(CC_data_iterated);
                    end
            end
        end
        
        function [RDI, CC_ws, CC_bs, RDI_control] = computeStability(obj, RespMat)
            

            if ~strcmp(obj.subsample_flag, 'Yes')
                RDI = zeros(obj.num_sessions-1, obj.num_cells);   %Instability Indicies, 1 is D7, 2 is D14, etc
                CC_ws = zeros(obj.num_sessions-1, obj.num_cells);
                CC_bs = zeros(obj.num_sessions-1, obj.num_cells);
                RDI_control = zeros(1, obj.num_cells);       % compare halves of trials in first session for stability control curve

                for ss = 2:obj.num_sessions
                    fprintf('Calculating Stability, session %d of %d\n', ss, obj.num_sessions);
                    [RDI(ss-1, :), CC_ws(ss-1, :), CC_bs(ss-1, :)] = ...
                        obj.stabilityBtwSessions(RespMat(:, :, :, 1), RespMat(:, :, :, ss));           %turn this into a method instead of a function
                end

                for qq = 1:10       %iterate 10 times and take the average
                    fprintf('Calculating chance stability...\n');
                    trialshuffle = randperm(size(RespMat, 1));
                    trialselect1 = trialshuffle(1:floor(size(RespMat, 1)/2));
                    trialselect2 = trialshuffle(floor(size(RespMat, 1)/2)+1:end);

                    half1 = RespMat(trialselect1, :, :, 1);
                    half2 = RespMat(trialselect2, :, :, 1);
                    [temp(qq, :), ~, ~] = obj.stabilityBtwSessions(half1, half2);
                end
                RDI_control = mean(temp, 1);
            else
                iterations = 10;

                RDI_temp = zeros(obj.num_sessions-1, obj.num_cells, iterations);   %Instability Indicies, 1 is D7, 2 is D14, etc
                CC_ws_temp = zeros(obj.num_sessions-1, obj.num_cells, iterations);
                CC_bs_temp = zeros(obj.num_sessions-1, obj.num_cells, iterations);
                RDI_control_temp = zeros(obj.num_cells, iterations);       % compare halves of trials in first session for stability control curve
                
                for rr = 1:iterations
                    trialshuffle = randperm(size(RespMat, 1));
                    trialselect = trialshuffle(1:obj.subsample);
                    RespMat_subsampled = RespMat(trialselect, :, :, :);

                    for ss = 2:obj.num_sessions
                        fprintf('Calculating Stability, session %d of %d\n', ss, obj.num_sessions);
                        [RDI_temp(ss-1, :, rr), CC_ws_temp(ss-1, :, rr), CC_bs_temp(ss-1, :, rr)] = ...
                            obj.stabilityBtwSessions(RespMat_subsampled(:, :, :, 1), RespMat_subsampled(:, :, :, ss));           %turn this into a method instead of a function
                    end

                    for qq = 1:10       %iterate 10 times and take the average
                        fprintf('Calculating chance stability...\n');
                        trialshuffle = randperm(size(RespMat_subsampled, 1));
                        trialselect1 = trialshuffle(1:floor(size(RespMat_subsampled, 1)/2));
                        trialselect2 = trialshuffle(floor(size(RespMat_subsampled, 1)/2)+1:end);

                        half1 = RespMat_subsampled(trialselect1, :, :, 1);
                        half2 = RespMat_subsampled(trialselect2, :, :, 1);
                        [temp(qq, :), ~, ~] = obj.stabilityBtwSessions(half1, half2);
                    end
                    RDI_control_temp(:, rr) = mean(temp, 1);
                end
                RDI = nanmean(RDI_temp, 3);
                CC_ws = nanmean(CC_ws_temp, 3);
                CC_bs = nanmean(CC_bs_temp, 3);
                RDI_control = nanmean(RDI_control_temp, 2);
            end
        end
        
        function [RDI, CC_ws, CC_bs] = stabilityBtwSessions(obj, Resp_Session1, Resp_Session2)

            % Resp matrices are [repeats x frames x cells]

            numReps = size(Resp_Session1,1);
            CC_ws = zeros(1,obj.num_cells); 
            CC_bs = zeros(1,obj.num_cells);
            switch obj.CC_method
                case 'Average'
                    for n = 1:obj.num_cells
                        ws1_alltrials = zeros(1,numReps);
                        ws2_alltrials = zeros(1,numReps);
                        bs1_alltrials = zeros(1,numReps);
                        bs2_alltrials = zeros(1,numReps);
                        for t = 1:numReps
                            temp = corrcoef(Resp_Session1(t,:,n),nanmean(Resp_Session1(1:numReps~=t,:,n)));
                            ws1_alltrials(t) = temp(2);
                            temp = corrcoef(Resp_Session2(t,:,n),nanmean(Resp_Session2(1:numReps~=t,:,n)));
                            ws2_alltrials(t) = temp(2);
                            temp = corrcoef(Resp_Session1(t,:,n),nanmean(Resp_Session2(1:numReps~=t,:,n)));
                            bs1_alltrials(t) = temp(2);
                            temp = corrcoef(Resp_Session2(t,:,n),nanmean(Resp_Session1(1:numReps~=t,:,n)));
                            bs2_alltrials(t) = temp(2);
                        end
                        CC_ws_alltrials = mean([ws1_alltrials; ws2_alltrials]);
                        CC_bs_alltrials = mean([bs1_alltrials; bs2_alltrials]);
                        CC_ws(n) = max(nanmean(CC_ws_alltrials),0);
                        CC_bs(n) = max(nanmean(CC_bs_alltrials),0);
                    end
                case 'Evenodd'
                    for n = 1:obj.num_cells
                        Eventrace1 = nanmean(Resp_Session1(2:2:end, :, n), 1);
                        Oddtrace1 = nanmean(Resp_Session1(1:2:end, :, n), 1);
                        Eventrace2 = nanmean(Resp_Session2(2:2:end, :, n), 1);
                        Oddtrace2 = nanmean(Resp_Session2(1:2:end, :, n), 1);
                        ws_1 = corr(Eventrace1', Oddtrace1');
                        ws_2 = corr(Eventrace2', Oddtrace2');
                        bs_1 = corr(Eventrace1', Oddtrace2');
                        bs_2 = corr(Eventrace2', Oddtrace1');
                        CC_ws(n) = max(mean([ws_1; ws_2]), 0);
                        CC_bs(n) = max(mean([bs_1; bs_2]), 0);
                    end
                case 'Random'
                    if numReps > 8
                        iterations = 100;
                    else
                        iterations = 50;
                    end
                    for n = 1:obj.num_cells
                        CC_ws_iterated = zeros(1, iterations);
                        CC_bs_iterated = zeros(1, iterations);
                        for jj = 1:iterations
                            trialshuffle = randperm(numReps);
                            trialselect1 = trialshuffle(1:floor(numReps/2));            %random subset of trials 1
                            trialselect2 = trialshuffle(floor(numReps/2)+1:end);        %random subset of trials 2
                            ses1_randavg1 = mean(Resp_Session1(trialselect1, :, n), 1); %avg of random subset 1, ses 1
                            ses1_randavg2 = mean(Resp_Session1(trialselect2, :, n), 1); %avg of random subset 2, ses 1
                            ses2_randavg1 = mean(Resp_Session2(trialselect1, :, n), 1); %avg of random subset 1, ses 2
                            ses2_randavg2 = mean(Resp_Session2(trialselect2, :, n), 1); %avg of random subset 2, ses 2
                            ws_1 = corr(ses1_randavg1', ses1_randavg2');
                            ws_2 = corr(ses2_randavg1', ses2_randavg2');
                            bs_1 = corr(ses1_randavg1', ses2_randavg2');
                            bs_2 = corr(ses1_randavg2', ses2_randavg1');  
                            CC_ws_iterated(jj) = max(nanmean([ws_1; ws_2]), 0);
                            CC_bs_iterated(jj) = max(nanmean([bs_1; bs_2]), 0);
                        end
                        CC_ws(n) = nanmean(CC_ws_iterated);
                        CC_bs(n) = nanmean(CC_bs_iterated);
                    end
            end

            if strcmp(obj.rectify_flag, 'Yes')
                CC_ws(CC_ws <= 0) = obj.CC_FUDGE;
            end

            RDI = (CC_ws-CC_bs)./(CC_ws+CC_bs);
        end
        
        %Getter
        function StabilityData = getStabilityData(obj)
            StabilityData = obj.StabilityData;
        end
        
        %Setters
        function setStabilityData(obj, input)
            obj.StabilityData = input;
        end
        
        function setCCmethod(obj, input)
            obj.CC_method = input;
        end

        function setSubsampleflag(obj, input)
            obj.subsample_flag = input;
        end
    end
end
        