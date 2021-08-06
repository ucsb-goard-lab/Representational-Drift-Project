
function [RFmapping_results] = ReceptiveField_analysis(fn_alt, pn_alt, fn_azi, pn_azi, save_flag)
% further basic processing and characterization of receptive field mapping stimulus responses

%% Written KS
%% Updated TM 190718
if nargin == 0
    disp('Choose your horizontal (altitude) data...')
    [fn_alt,pn_alt] = uigetfile('.mat');
    disp('Choose your vertical (azimuth) data...')
    [fn_azi,pn_azi] = uigetfile('.mat');
    save_flag = questdlg('Save data?', 'Query','Yes', 'No', 'Yes');
end

alt_data = importdata(fn_alt);
azi_data = importdata(fn_azi);

altitude = alt_data.RespVec;
azimuth = azi_data.RespVec;

% isVisuallyResponsive = alt_data.isVisuallyResponsive | azi_data.isVisuallyResponsive;

% Initial pass, let's just mean everything and turn it into a heatmap...

altitude_m = squeeze(mean(mean(altitude,2),4));
azimuth_m = squeeze(mean(mean(azimuth,2),4));

for ii = 1:size(altitude_m,1)
    [azi_fit{ii},azigof] = fit([1:size(azimuth_m,2)]',azimuth_m(ii,:)','gauss1','Upper',[20,size(azimuth_m,2),20],'Lower',[0,1,0]);
    azi_pref(ii) = azi_fit{ii}.b1;
    azi_width(ii) = azi_fit{ii}.c1;
    azi_rsquare(ii) = azigof.rsquare;
    
    [alt_fit{ii},altgof] = fit([1:size(altitude_m,2)]',altitude_m(ii,:)','gauss1','Upper',[20,size(altitude_m,2),20],'Lower',[0,1,0]);
    alt_pref(ii) = alt_fit{ii}.b1;
    alt_width(ii) = alt_fit{ii}.c1;
    alt_rsquare(ii) = altgof.rsquare;
end

for ii = 1:size(altitude,1)
    alt_p(ii) = anova1(squeeze(mean(altitude(ii,:,:,:),2))',[],'off');
    azi_p(ii) = anova1(squeeze(mean(azimuth(ii,:,:,:),2))',[],'off');
end


%edge discarding

% isEdgeAzi = false(1,size(azimuth_m,1));
% isEdgeAlt = false(1,size(altitude_m,1));
%     for ii = 1:size(azimuth_m,1)
%         if round(azi_fit{ii}.b1,3) == 1 || round(azi_fit{ii}.b1,3) == size(azimuth_m,2)
%         isEdgeAzi(ii) = true;
%         end
%         
%         if round(alt_fit{ii}.b1,3) == 1 || round(alt_fit{ii}.b1,3) == size(altitude_m,2)
%         isEdgeAlt(ii) = true;
%         end
%     end
isTooNarrow = (azi_width<1) | (alt_width<1);
% isEdge = isEdgeAzi | isEdgeAlt;
% isTun = (alt_p<0.05) | (azi_p<0.05);

isFit = (alt_rsquare>0.5) | (azi_rsquare>0.5);

isSpatiallyTuned = isFit & ~isTooNarrow;
isSpatiallyTunedAlt = alt_rsquare>0.5 & ~alt_width<2;
isSpatiallyTunedAzi = azi_rsquare>0.5 & ~azi_width<2;

disp(['Fraction spatially tuned: ' num2str(mean(isSpatiallyTuned))])


%% extracting centroids
%load your data file...
numCells = size(azi_data.cellMasks,2);
for c = 1:numCells
    currCellMask = azi_data.cellMasks{c};
   temp = regionprops(poly2mask(currCellMask(:,1),currCellMask(:,2),760,760),'centroid');
   roi_centroids(c,:) = temp.Centroid;
end
   
isANOVA = alt_p<0.01 | azi_p<0.01;

isSpat = isSpatiallyTuned & isANOVA;
x_locations = roi_centroids(:,1);
y_locations = roi_centroids(:,2);
try
    alt_corr = corr(alt_pref(isSpat)',x_locations(isSpat));
    azi_corr = corr(azi_pref(isSpat)',y_locations(isSpat));
catch
    alt_corr = NaN;
    azi_corr = NaN;
end

if strcmp(save_flag, 'Yes')
	save RFmapping_results.mat azi_pref alt_pref altitude azimuth isSpatiallyTuned alt_p azi_p azi_fit alt_fit roi_centroids alt_corr azi_corr
end

RFmapping_results.alt_pref = alt_pref;
RFmapping_results.azi_pref = azi_pref;
RFmapping_results.altitude = altitude;
RFmapping_results.azimuth = azimuth;
RFmapping_results.isSpatiallyTuned = isSpatiallyTuned;
RFmapping_results.isSpatiallyTunedAlt = isSpatiallyTunedAlt;
RFmapping_results.isSpatiallyTunedAzi = isSpatiallyTunedAzi;
RFmapping_results.alt_p = alt_p;
RFmapping_results.azi_p = azi_p;
RFmapping_results.alt_fit = alt_fit;
RFmapping_results.azi_fit = azi_fit;
RFmapping_results.roi_centroids = roi_centroids;
RFmapping_results.alt_corr = alt_corr;
RFmapping_results.azi_corr = azi_corr;


end