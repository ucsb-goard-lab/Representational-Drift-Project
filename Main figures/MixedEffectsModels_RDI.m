function [pVal, F, DF1, DF2, RDI_stim1_mean, RDI_stim1_SE, RDI_stim2_mean, RDI_stim2_SE] = MixedEffectsModels_RDI()
%% Load data
RDIdata = importdata('RDIdata.mat');
% numStim = 2;                                   % number of stimuli 
numMice = length(RDIdata);                     % calculate number of mice
%% Correct for non-consecutive weeks
fields = fieldnames(RDIdata(1));
stim_1 = fields{1};
stim_2 = fields{2};
for m = 1:numMice
    weeks = RDIdata(m).weeks;
    curr_numSess = size(RDIdata(m).present1, 1);
    if curr_numSess < find(weeks, 1, 'last')     % if there's a mismatch between number of sessions and the latest week in the real-time weeks vector
        RDIdata(m).(stim_1) = alignWeeks(RDIdata(m).(stim_1), weeks);     % insert dummy data to shift data such that it lines up with real-time weeks
        RDIdata(m).(stim_2) = alignWeeks(RDIdata(m).(stim_2), weeks);
        RDIdata(m).present1 = alignWeeks(RDIdata(m).present1, weeks);
        RDIdata(m).present2 = alignWeeks(RDIdata(m).present2, weeks);
    end
end

for m = 1:numMice           
    numSess(m) = size(RDIdata(m).present1,1);   % calculate sessions per mouse
%     numSamp(m) = sum(sum(RDIdata(m).present)); % calculate samples per mouse
end
pVal = zeros(1, max(numSess));
F = zeros(1, max(numSess));
%% Collate data
for s = 1:max(numSess)
    
    % Initialize matrices
    RDI = [];             % Representational Drift Index (Response variable)
    StimID = [];          % Stimulus: 0 = first stimulus, 1 = second stimulus (Fixed effect)
    SessionID = [];       % Session number (Fixed Effect)
    MouseID = [];         % Mouse ID (Random effect)
    
    for m = 1:numMice
        if RDIdata(m).weeks(s) == 1      % if this mouse has data for this (real-time) week
            % Determine which cells are present for this session, separately for each "stimulus" input in case we have a different number of neurons for each
            % present1 and present2 will be identical if using the same neurons
            present_idx_1 = RDIdata(m).present1(s,:);
            present_idx_2 = RDIdata(m).present2(s,:);
            % Add RDIs for both conditions for present cells
            RDI_stim1 = RDIdata(m).(stim_1)(s,logical(present_idx_1));
            RDI_stim2 = RDIdata(m).(stim_2)(s,logical(present_idx_2));
            RDI = [RDI; RDI_stim1'; RDI_stim2'];
            % Add Stim IDs for present cells
            StimID = [StimID; zeros(sum(present_idx_1),1); ones(sum(present_idx_2),1)];
            % Add Mouse IDs fro present cells
            MouseID = [MouseID; m*ones(sum(present_idx_1)+sum(present_idx_2),1)];

            % Extract RDI data by mouse for plotting
            RDI_stim1_MouseAvg(m) = nanmean(RDI_stim1);
            RDI_stim2_MouseAvg(m) = nanmean(RDI_stim2);
        end
    end
    
    % Calculate mean RDI and standard error across mice (for plotting)
    RDI_stim1_mean(s) = nanmean(RDI_stim1_MouseAvg);
    RDI_stim1_SE(s) = nanstd(RDI_stim1_MouseAvg)/sqrt(length(unique(MouseID)));
    RDI_stim2_mean(s) = nanmean(RDI_stim2_MouseAvg);
    RDI_stim2_SE(s) = nanstd(RDI_stim2_MouseAvg)/sqrt(length(unique(MouseID)));

    
    % Put data into table
    tbl = table(RDI,StimID,MouseID,'VariableNames',{'RDI','StimID','MouseID'});
    
    % Fit linear mixed-effects model
    lme = fitlme(tbl,'RDI~1+StimID+(StimID|MouseID)','FitMethod','REML');
    
    % Show all stats for each session
    disp(['Session ' num2str(s) '/' num2str(max(numSess))])
%     pause
    
    % Extract p values
    [pVal(s),F(s),DF1(s),DF2(s)] = coefTest(lme);
    
end

%% Plot

% data
% figure
if length(RDI_stim1_mean) > 1
    errorbar([2:7],RDI_stim1_mean,RDI_stim1_SE,'color',[0 0.8 0.5],'linewidth',2)
    hold on
    errorbar([2:7],RDI_stim2_mean,RDI_stim2_SE,'color',[0.8 0 0.5],'linewidth',2)
%     plot([1:8],zeros(1,8),'k:')
else
    errorbar(2,RDI_stim1_mean,RDI_stim1_SE,'color',[0 0.8 0.5],'linewidth',2)
    hold on
    errorbar(2,RDI_stim2_mean,RDI_stim2_SE,'color',[0.8 0 0.5],'linewidth',2)
%     plot([1:3],zeros(1,3),'k:')
end
axis square

% asterisks
for s = 1:max(numSess)
    if pVal(s)<0.05
        plot(s+1,max(RDI_stim1_mean(s)+RDI_stim1_SE(s),RDI_stim2_mean(s)+RDI_stim2_SE(s))+0.05,'k*')
    end
end

% make pretty
axis([1.5 7.5 -0.1 0.5])
set(gcf,'color',[1 1 1])

end

function corrected = alignWeeks(indata, recording_weeks)
% this subfunction assumes indata with the format: [sessions x cells]
num_cells = size(indata, 2);
insert_vec = NaN(1, num_cells);             % insert a vector of NaNs in the shift index location, and move everything one over

shift_idx = find(diff(recording_weeks) == -1) + 1;       % find all locations where data needs to be shifted
num_shifts = length(shift_idx);

corrected = indata;
for tt = 1:num_shifts
    if shift_idx(tt) <= find(recording_weeks, 1, 'last')
        corrected = cat(1,  corrected(1:shift_idx(tt)-1, :), insert_vec, corrected(shift_idx(tt):end, :));
    end
end
    

end
