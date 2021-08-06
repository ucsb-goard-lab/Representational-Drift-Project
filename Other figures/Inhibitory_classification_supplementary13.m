%% Inhibitory neuron categorization


an = StabilityAnalyzer();
an.importData();
an.setQualityThreshold(3);
an.cellSelection;
sessions = an.getUse_sessions;
num_cells = an.num_cells;
num_sessions = an.num_sessions;

sessions(isnan(sessions)) = 1;            % nan values converted to zeros for logical purposes
present = sum(sessions, 2) == an.num_sessions;

oridata = an.getOrientationData;

thresh = 0.4;
osi = oridata.osiMat;
avgosi = zeros(1, num_cells);
stdosi = zeros(1, num_cells);
for ii = 1:size(osi, 2)
    avgosi(ii) = mean(osi(logical(sessions(ii, :)), ii));
    stdosi(ii) = std(osi(logical(sessions(ii, :)), ii));
end
sharp = avgosi >= thresh;
broad = avgosi < thresh;



% running prep for nested stats to plug into LMEM
RDIdata = struct();
mouse = 0;


stim = 'PDG';
qual = 3; 

% navigate to each directory and run the following (with appropriate 'weeks' indicator for the given mouse)
weeks = [1 1 1 1 1 1];          % week indicator, starting on week 2. e.g. [1 1 0 1 1 1] means data (subsequent to first session) comes from weeks 2, 3, 5, 6, 7. 
mouse = mouse + 1;
load('StabilityData.mat');
load('RoiINFO.mat');
sharp_dual = RoiINFO.PDG_Responsive_thresh & RoiINFO.NatMov_Responsive_thresh & RoiINFO.quality >= qual & sharp;
broad_dual = RoiINFO.PDG_Responsive_thresh & RoiINFO.NatMov_Responsive_thresh & RoiINFO.quality >= qual & broad;
presence = RoiINFO.presence;

stim_sharp = StabilityData.(stim).RDI(:, sharp_dual);
stim_broad = StabilityData.(stim).RDI(:, broad_dual);
present_sharp = presence(sharp_dual, 2:end)';
present_broad = presence(broad_dual, 2:end)';

RDIdata(mouse).stim_sharp = stim_sharp;
RDIdata(mouse).stim_broad = stim_broad;
RDIdata(mouse).present1 = present_sharp;             % one presence indicator for each stimuli in case we are using different neurons and thus have different sample sizes between them
RDIdata(mouse).present2 = present_broad;
RDIdata(mouse).weeks = weeks;

mouse



