% Master script for further processing DF/F data into format for analysis (format in which the data for this study is provided). 
% Prodcues data organized into different structures for ease of use.
%
%   NOTE: "MOV" == "NatMov"
%
% RespData = [repeats, frames, cells, sessions], contains response matrices for both stimuli for all neurons across all sessions 
%       - RespMat_Full = full response matrix
%       - RespMat_onTime = response matrix excluding pre-stimulus gray screen period (MOV only)
%       
% StabilityData: contains reliability and stability information (RDI & constituent components) for every ROI
%       - CCs = [sessions x cells] reliability values for every ROI across sessions
%       - RDI = [sessions x cells] RDI values for every ROI (from 2nd session to last session)
%       - CC_ws = [sessions x cells] within session CC value for RDI calculation
%       - CC_bs = [sessions x cells] between session CC value for RDI calculation
%       - RDI_control = [1 x cells] control RDI value for every ROI (RDI calculated using session 1 data)
%
% RoiINFO: contains ROI visual responsiveness indicators for both stimuli, as well as a quality rating (and by-session presence indicator) determined through visual inspection
%       - quality = [1 x cells] ROI quality rating determined by visual inspection
%       - presence = [cells x sessions] logical indicator for neuron presence on a given session
%       - PDG_isTuned = [1 x cells] logical indicator for neurons orientation tuned on all sessions
%       - PDG_Responsive_thresh = [1 x cells] logical indicator for neurons responsive to the PDG stimulus
%       - NatMov_Responsive_thresh = [1 x cells] logical indicator for neurons responsive to the MOV stimulus
%       - other fields are not relevant
%
% OrientationData: contains orientation tuning information of ROIs across sessions for the PDG stimulus
%       - isTuned = [sessions x cells] logical indicaotr for orientation tuned neurons
%       - osiMat = [sessions x cells] OSI values for all ROIs across sessions
%       - zScoreMat = [sessions x cells] zScore of average onset responses to preferred stimuli
%       - oriResp = [sessions x cells x orientations] average responses to each orientation across sessions
%       - oriPref = [sessions x cells] preferred orientation across sessions
%       - osiAvg = [sessions x 1] average OSI 
%       - zScoreAvg = [sessions x 1] average zScore of onset responses to preferred stimuli
%       - osiSEM & zScoreSEM = [sessions x 1] corresponding standard errors for osiAvg and zScoreAvg
%

% initialize mouse object
tdm091 = ChronicImaging_Mouse;

% set method for computing reliability ('Random' = CC between trial averaged activity from two random halves of trials, iterated and averaged)
CC_method = 'Random';

% create and run processor objects for each stimulus
pdg = ProcessorController('PDG');
pdg.processor.setCCmethod(CC_method);
pdg.processor.setSubsampleflag('No');           % no subsampling for PDG
pdg.processor.setOriflag('Yes');                % YES to perform orientation tuning analysis
pdg.processor.setDataType('DFF');               % is input data DFF or spikes?
pdg.runProcessor();

natmov = ProcessorController('NatMov');
natmov.processor.setCCmethod(CC_method);
natmov.processor.setSubsampleflag('Yes');            %YES for DFF data, NO for spike data (don't need any CC data or stability data for spike data, only need the Resp matrix for later event detection)
natmov.processor.setDataType('DFF');        
natmov.runProcessor(); 

% PDG stimulus with no interstimulus periods
pdgcont = ProcessorController('PDGcont');
pdgcont.processor.setCCmethod(CC_method);
pdgcont.processor.setSubsampleflag('Yes');             
pdgcont.processor.setDataType('DFF');
pdgcont.runProcessor();

% NatMov stimulus with interstimulus periods
natmovdisc = ProcessorController('NatMovdisc');
natmovdisc.processor.setCCmethod(CC_method);
natmovdisc.processor.setSubsampleflag('No');
natmovdisc.processor.setDataType('DFF');
natmovdisc.processor.setOriflag('No');
natmovdisc.runProcessor();



% import processor objects into mouse object
tdm091.importObjects('PDG', pdg, 'NatMov', natmov);
% tdm091.importObjects('PDG', pdg, 'NatMov', natmov, 'PDGcont', pdgcont, 'NatMovdisc', natmovdisc);
% tdm091.importObjects('PDGcont', pdgcont, 'NatMovdisc', natmovdisc);

% or import previously saved data structures 
tdm091.importData();

% Evaluate ROIs for quality and visual responsiveness
%   - finds responsive cells (select 'PDG_Responsive_thresh', 'NatMov_Responsive_thresh', 'isTuned')
%   - performs visual ROI inspection (select 'Visual_Inspection_A' for inspection program, and then 'Visual_Inspection_upload' when finished to upload saved cellInfo structure to mouse object)
tdm091.evaluateROIs;

% save desired data
tdm091.saveData('all');


% plot example RDI curves (average +- sem across neurons, note: RDI curves in actual figures are produced by MixedEffectsModels_RDI function)
test = StabilityAnalyzer;
test.importData();
% test.addData();               % add additional mice one at a time for pooling data
test.setQualityThreshold(3);    % set quality threshold
test.cellSelection;             % select inclusion criteria (e.g. 'quality', 'presence', 'PDG_Responive_thresh', 'NatMov_Responsive_thresh' for quality neurons responsive to both stimuli)

figure;
test.RDIplotter('PDG', 'Average');       % select inclusion criteria for the PDG curve
hold on
test.RDIplotter('NatMov', 'Average');    % select inclusion criteria for the MOV curve
% tdm092a.RDIplotter('MultiMov', 'Average');
% test.RDIplotter('PDGcont', 'Average');
% test.RDIplotter('NatMovdisc', 'Average');
legend off
set(gcf, 'Position', [800 300 900 900]);



% receptive field stuff, same process as above
tdm091 = ChronicImaging_Mouse;
rf = ProcessorController('RF');
rf.runProcessor();
tdm091.importObjects('RF', rf);
tdm091.saveData('RFdata');

% visualizing responses, same process as above. Import first field and add additional fields if desired
tdm0091v = ResponseVisualizer();
tdm0091v.importData();
tdm0091v.addData();
tdm0091v.setQualityThreshold(3);
tdm0091v.cellSelection;

tdm0091v.plotPDGvNatMov_traces2(1, 'minimum');          % first arg: starting cell, second arg: normalization type
tdm0091v.plotTempoStructure_traces(1, 'minimum');
tdm0091v.plotPDGvNatMov_avgs(1);
