% Prepares data for use in MixedEffectsModels_RDI function for generating RDI curves

%% for normal PDG vs MOV comparisons
% initialize RDIdata and mouse once
RDIdata = struct();
mouse = 0;

% Select stimuli to compare. LMEM will be reading different labels than the fieldnames in the original data.
stim1_in = 'PDG';            % readable from data files. Change these based on what you want to compare using stimulus names in original datafiles
stim2_in = 'NatMov';
stim1_out = 'PDG';              % readable by the LME tester. These labels are arbitrary and only for LMEM to read. Keep them the same and know what you are comparing
stim2_out = 'MOV';
qual = 3; 

% navigate to each directory and run the following (with appropriate 'weeks' indicator for the given mouse)
weeks = [1 1 1 1 1 1];          % week indicator, starting on week 2. e.g. [1 1 0 1 1 1] means data (subsequent to first session) comes from weeks 2, 3, 5, 6, 7. 
mouse = mouse + 1;
load('StabilityData.mat');
load('RoiINFO.mat');
dualresp = RoiINFO.PDG_Responsive_thresh & RoiINFO.NatMov_Responsive_thresh & RoiINFO.quality >= qual;
movresp = RoiINFO.NatMov_Responsive_thresh & RoiINFO.quality >= qual;
pdgresp = RoiINFO.PDG_Responsive_thresh & RoiINFO.quality >= qual;
all = RoiINFO.quality >= qual;
presence = RoiINFO.presence;

stim1 = StabilityData.(stim1_in).RDI(:, dualresp);
stim2 = StabilityData.(stim2_in).RDI(:, dualresp);
present1 = presence(dualresp, 2:end)';
present2 = presence(dualresp, 2:end)';

RDIdata(mouse).(stim1_out) = stim1;
RDIdata(mouse).(stim2_out) = stim2;
RDIdata(mouse).present1 = present1;             % one presence indicator for each stimuli in case we are using different neurons and thus have different sample sizes between them
RDIdata(mouse).present2 = present2;
RDIdata(mouse).weeks = weeks;

mouse

%% multimov comparisons
RDIdata = struct();
mouse = 0;
stim1_in = 'PDG';            % readable from data files
stim2_in = 3;                 % which scramble to use, 1 = 0%, 2 = 50%, 3 = 100%
stim1_out = 'PDG';              % readable by the LME tester
stim2_out = 'MOV';

weeks = [1 1 1 1 1 1];
mouse = mouse + 1;
load('StabilityData.mat');
load('RoiINFO.mat');
dualresp = RoiINFO.PDG_Responsive_thresh & RoiINFO.MultiMov_Responsive_thresh & RoiINFO.quality >= 3;
movresp = RoiINFO.MultiMov_Responsive_thresh & RoiINFO.quality >= 3;
pdgresp = RoiINFO.PDG_Responsive_thresh & RoiINFO.quality >= 3;
all = RoiINFO.quality >= 3;
presence = RoiINFO.presence;

stim1 = StabilityData.(stim1_in).RDI(:, dualresp);
stim2 = StabilityData.MultiMov.RDI{stim2_in}(:, dualresp);
present1 = presence(dualresp, 2:end)';
present2 = presence(dualresp, 2:end)';

RDIdata(mouse).(stim1_out) = stim1;
RDIdata(mouse).(stim2_out) = stim2;
RDIdata(mouse).present1 = present1;             % one presence indicator for each stimuli in case we are using different neurons and thus have different sample sizes between them
RDIdata(mouse).present2 = present2;
RDIdata(mouse).weeks = weeks;

mouse