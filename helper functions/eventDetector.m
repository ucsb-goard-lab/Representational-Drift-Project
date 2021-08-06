function [eventvectors] = eventDetector(Resp, F0, Resp_sp, F0_sp, min_length, sig_thresh, session_thresh)
    % Inputs:
    % Resp = [repeats x frames x cells x sessions] 
    % F0 = [repeats x cells x sessions] average baseline fluorescence value to be used to normalize Resp
    %  _sp modifier is the above for spiking data: deconvolved spike data is included for refinement of event periods determined by calcium data
    % min_length = minimum length of event periods in frames
    % sig_thresh = significance threshold (alpha value) for stat test to determine visually responsive frames
    % session_thresh = minimum number of sessions for which a frame must pass visual responsiveness test to be considered event frames
    % 
    % Outputs:
    % eventvectors = [cells x frames] vectors describing event periods for each neuron
    %
    % Written 200730 TM

    % test for events on each session separately, if an event is detected on any session it will be counted for all sessions
    
    % initialize
    eventvectors = zeros(size(Resp, 3), size(Resp, 2));
    
    spike_cleaning = 1;     % use spike data to clean up even periods
    gap_thresh = 3;         % permissible gap length in between events (frames)

    % Normalize trials of each session with average baseline fluorescence on respective trial and session
    normResp = zeros(size(Resp));       % for DFF data
    for ii = 1:size(Resp, 3)
        for kk = 1:size(Resp, 4)
            normResp(:, :, ii, kk) = Resp(:, :, ii, kk) - F0(:, ii, kk);
        end
    end

    % Detect events         
    for ii = 1:size(normResp, 3)
        ii

        % Detect significant activity
        curr_vector = zeros(1, size(normResp, 2));          % [1 x frames] temp logical vector describing significant frames for cell ii
        for ff = 1:size(normResp, 2)
            currdetection = zeros(1, size(normResp, 4));    % [1 x sessions] temp logical vector describing significance of frame ff on each session
            for kk = 1:size(Resp, 4)
                curr_cell = squeeze(normResp(:, ff, ii, kk));
                p = signrank(curr_cell, 0, 'tail', 'right');
                if p < sig_thresh
                    currdetection(kk) = 1;              
                end
            end
            if sum(currdetection) >= session_thresh     % if frame is found to be significant for at least session_thresh sessions, indicate significance in frame vector
                curr_vector(ff) = 1;
            end
        end

        % Clean up: delete isolated events
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
            if ele_counter < min_length                         % number of times we incremented our contiguity counter corresponds to the length of the 'event' we're looking at
                curr_vector(curr_idx+1:curr_idx+ele_counter) = 0;       % if it's less than a minimum length, set all the values of this 'event' to 0
            end
        end

        % Clean up: combine very close events
        switches = diff(curr_vector);
        switch_idcs = find(switches == -1);
        for qq = 1:length(switch_idcs)
            curr_idx = switch_idcs(qq);
            ele_counter = 1;
            end_flag = 0;
            while ~end_flag
                if (curr_idx + ele_counter) <= length(switches)
                    curr_val = switches(curr_idx + ele_counter);

                    if curr_val == 0
                        ele_counter = ele_counter + 1;
                    else
                        end_flag = 1;
                    end
                else
                    end_flag = 1;
                end
            end
            if ele_counter < gap_thresh
                curr_vector(curr_idx+1:curr_idx+ele_counter) = 1;
            end
        end

        eventvectors(ii, :) = curr_vector;
    end
  
    if spike_cleaning == 1
        % find periods where spike rate is greater than 2 sp/s
        catSpikes = [];
        for kk = 1:size(Resp_sp, 4)     % create PSTH
            catSpikes = cat(1, catSpikes, Resp_sp(:, :, :, kk));
        end

        psth = squeeze(mean(catSpikes, 1));
        for ii = 1:size(psth, 2)
            super(ii, :) = psth(:, ii) > 2;             % spikes per second
        end

        % for every cell, iterate through every event - if the end of the event is more than n (10) frames away from the last significant frame before it in super, set end of event through start of next super period to 0
        grace_period = 10;
        for ii = 1:size(eventvectors, 1)
            switches_events = diff(eventvectors(ii, :));            % -1s indicate the end frames of every event
            switches_super = diff(super(ii, :));                    % -1s indicate the end frame of every super spike threshold period

            event_ends = find(switches_events == -1);
            super_ends = find(switches_super == -1);
            super_starts = find(switches_super == 1) + 1;
            
            %some cleanup
            if length(super_starts) > length(super_ends)       % case in which event starts on first frame
                super_ends(end+1) = size(eventvectors, 2);
            elseif length(super_ends) > length(super_starts)       % case in which event stops on last frame
                super_starts = [1 super_starts];
            elseif length(super_starts) > 1                   % case in which both of the above happen in the same neuron
                if super_starts(1) > super_ends(1)
                    super_ends(end+1) = size(eventvectors, 2);
                    super_starts = [1 super_starts];
                end
            end
            
            for ee = 1:length(super_ends)
                curr_super_end = super_ends(ee);
                try
                    next_super_start = super_starts(ee+1);    
                catch
                    next_super_start = size(eventvectors, 2);
                end
                following_event_ends = event_ends(event_ends > curr_super_end);
                try                                     % if there are no more event ends
                    next_event_end = following_event_ends(1);
                catch
                    next_event_end = size(eventvectors, 2);
                end
                
                if next_event_end - curr_super_end > grace_period && next_super_start - curr_super_end > grace_period
                    eventvectors(ii, curr_super_end+grace_period:next_super_start-1) = 0;
                end

            end
        end 
    end
    
    % one more round of cleanup (turn this into a subfunction maybe)
    for ii = 1:size(normResp, 3)
        curr_vector = eventvectors(ii, :);
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
            if ele_counter < min_length                         % number of times we incremented our contiguity counter corresponds to the length of the 'event' we're looking at
                curr_vector(curr_idx+1:curr_idx+ele_counter) = 0;       % if it's less than a minimum length, set all the values of this 'event' to 0
            end
        end
        eventvectors(ii, :) = curr_vector;
    end

end
