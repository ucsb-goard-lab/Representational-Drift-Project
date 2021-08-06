
function [data, kernel] = spikeInference(filename,save_flag)

% Infers spike rates from calcium imaging data
% Expects input file with a fieldname 'DFF', representing calcium data in the form [neurons x frames]

%% Deconvolution code citations:

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

if nargin == 0
    [filename,pathname] = uigetfile('.mat');
    load(filename);
    save_flag = questdlg('Save file?','Dialog box','Yes','No','Yes');
elseif nargin == 1
    save_flag = questdlg('Save file?','Dialog box','Yes','No','Yes');
end

load(filename);

dffDeconv = zeros(size(data.DFF));

for n = 1:size(data.DFF, 1)
subroutine_progressbar(n/size(data.DFF,1));
    % get trace and run deconvolution
    n
    trace = data.DFF(n,:);
    [denoised,spikes,opt,kernel] = deconvolveCa(trace, 'ar1' ,'foopsi', 'optimize_pars');
    dffDeconv(n,:) = spikes;
    
end
subroutine_progressbar(1);close all;

data.spikes = dffDeconv;


if strcmp(save_flag,'Yes')
save(filename,'data');
end
end

