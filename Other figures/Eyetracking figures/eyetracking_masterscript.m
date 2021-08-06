% initialize the processor
eye = Eye_Processor;
% first, imports data, crops and resizes each session's video to common coordinates based on user input, and then detects pupil
eye.runA('Yes');
% run cropper or detector separately if necessary
eye.runCropper();
eye.runPupilDetector('No');
% second, packages data for analysis 
eye.runB();
eye.saveData();

post66 = EyeAnalyzer;
post66.importData();
post66.drawCentroidMap();




%% eyemovement figure
% start with first field
week_vec = [1 1 1 1 1 1 1];         % (for example)
[WS_cat, WS_avg] = eyemovement_comparisons(week_vec, 'r');      % week_vec = data indicator (e.g. [1 1 1 0 1 1 1] indicates data corresponds to weeks 1, 2, 3, 5, 6, 7)

% call function again with subsequent fields, input previously generated data structs
week_vec = [1 1 1 0 1 1 1];         % (for example)
[WS_cat, WS_avg] = eyemovement_comparisons(week_vec, 'b', WS_cat, WS_avg);          % use a different color for each field


%% pupil size figure
% start with first field
week_vec = [1 1 1 1 1 1 1];
[A_cat, A_avg, trial_A_cat, gains_cat] = pupilsize_comparisons(week_vec, 'r');

% call function again with subsequent fields, input previously generated data structs
week_vec = [1 1 0 1 1 1];       % (for example)
[A_cat, A_avg, trial_A_cat, gains_cat] = pupilsize_comparisons(week_vec, 'b', A_cat, A_avg, trial_A_cat, gains_cat);



