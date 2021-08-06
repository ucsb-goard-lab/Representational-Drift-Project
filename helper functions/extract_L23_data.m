% extracts L2/3 neuron data from mice with multiple layers 
load('LayerVec.mat')
L23 = strcmp(LayerVec, 'L2/3');
L4 = strcmp(LayerVec, 'L4');
L5 = strcmp(LayerVec, 'L5');

load('OrientationData.mat')
OrientationData.isTuned(:, ~L23) = [];
OrientationData.osiMat(:, ~L23) = [];
OrientationData.zScoreMat(:, ~L23) = [];
OrientationData.oriResp(:, ~L23, :) = [];
OrientationData.oriPref(:, ~L23) = [];
OrientationData.osiAvg = mean(OrientationData.osiMat, 2);
OrientationData.zScoreAvg = mean(OrientationData.zScoreMat, 2);
OrientationData.osiSEM = std(OrientationData.osiMat, [], 2);
OrientationData.zScoreSEM = std(OrientationData.zScoreMat, [], 2);

load('RespData.mat')
RespData.PDG.RespMat_Full(:, :, ~L23, :) = [];
RespData.NatMov.RespMat_Full(:, :, ~L23, :) = [];
RespData.NatMov.RespMat_onTime(:, :, ~L23, :) = [];

load('RoiINFO.mat')
RoiINFO.quality(~L23) = [];
RoiINFO.presence(~L23, :) = [];
RoiINFO.all(~L23) = [];
RoiINFO.PDG_Responsive_thresh(~L23) = [];
RoiINFO.PDG_Responsive_shuffle(~L23) = [];
RoiINFO.NatMov_Responsive_thresh(~L23) = [];
RoiINFO.NatMov_Responsive_shuffle(~L23) = [];
RoiINFO.PDG_isTuned(~L23) = [];

load('StabilityData.mat')
StabilityData.PDG.CCs(:, ~L23) = [];
StabilityData.PDG.RDI(:, ~L23) = [];
StabilityData.PDG.CC_ws(:, ~L23) = [];
StabilityData.PDG.CC_bs(:, ~L23) = [];
StabilityData.PDG.RDI_control(:, ~L23) = [];
StabilityData.NatMov.CCs(:, ~L23) = [];
StabilityData.NatMov.RDI(:, ~L23) = [];
StabilityData.NatMov.CC_ws(:, ~L23) = [];
StabilityData.NatMov.CC_bs(:, ~L23) = [];
StabilityData.NatMov.RDI_control(~L23, :) = [];

save OrientationData OrientationData
save RespData RespData
save RoiINFO RoiINFO
save StabilityData StabilityData