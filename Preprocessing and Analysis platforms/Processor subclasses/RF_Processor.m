classdef RF_Processor < ResponseProcessor
    
    properties (Access = protected)
        RFdata              % data structure for receptive field stimulus experiments 
        pn_data             %pathname for rf data
        pn_stimdata         %pathname for rf stimulus data files
    end
    
    methods
        function obj = RF_Processor()
            obj.RFdata = struct();
        end
        
        function run(obj)
            obj.getRFfiles;
            obj.Stimdata = obj.getStimData('RF', obj.ref_session);
            
            [RespMat_alt, RespMat_azi] = obj.processorA('Yes');
            obj.RespData = structPacker([], RespMat_alt, 'altitude', RespMat_azi, 'azimuth');
            
            obj.RFdata = obj.processorB;
%             obj.RFdata.RespData = obj.RespData;
            
        end  
        
        function getRFfiles(obj)
            disp('Select all ALTITUDE DATA files.')
            [obj.filelist.altdata, obj.pn_data] = uigetfile('.mat', 'MultiSelect', 'on');
            disp('Select all AZIMUTH DATA files.')
            [obj.filelist.azidata, ~] = uigetfile('.mat', 'MultiSelect', 'on');
            disp('Select all ALTITUDE STIMULUS files.')
            [obj.filelist.altstimdata, obj.pn_stimdata] = uigetfile('.mat', 'MultiSelect', 'on');
            disp('Select all AZIMUTH STIMULUS files.')
            [obj.filelist.azistimdata, ~] = uigetfile('.mat', 'MultiSelect', 'on');
            cd(obj.pn_data)
            
            lengthList = length(obj.filelist.altdata);
            
            if isempty(obj.reference_number)
                obj.reference_number = 1;              %reference session defaults to 1
            end
            obj.num_sessions = lengthList;
            obj.ref_session = importdata(obj.filelist.altdata{obj.reference_number});
            % azi_refSession = importdata(obj.filelist.RF.azidata{obj.reference});
            obj.num_cells = length(obj.ref_session.cellMasks);
        end
        
        function [RespMat_alt, RespMat_azi] = processorA(obj, overwrite_flag)
            isVisuallyResponsive_alt = zeros(obj.num_cells, obj.num_sessions);
            isVisuallyResponsive_azi = zeros(obj.num_cells, obj.num_sessions);
            RespMat_alt = zeros(obj.num_cells, obj.Stimdata.numReps,  obj.Stimdata.altlocs,  obj.Stimdata.onTime*obj.Stimdata.framerate, obj.num_sessions);
            RespMat_azi = zeros(obj.num_cells, obj.Stimdata.numReps,  obj.Stimdata.azilocs,  obj.Stimdata.onTime*obj.Stimdata.framerate, obj.num_sessions);
            
            for kk = 1:obj.num_sessions
                fprintf('Pre-processing altitude file %d out of %d\n', kk, obj.num_sessions);
                ReceptiveField_preProcessing(obj.filelist.altstimdata{kk}, obj.pn_stimdata, obj.filelist.altdata{kk}, obj.pn_data, overwrite_flag);
                currdata = importdata(obj.filelist.altdata{kk});
                isVisuallyResponsive_alt(:, kk) = currdata.isVisuallyResponsive;
                RespMat_alt(:, :, :, :, kk) = currdata.RespVec;
                fprintf('Pre-processing azimuth file %d out of %d\n', kk, obj.num_sessions);
                ReceptiveField_preProcessing(obj.filelist.azistimdata{kk}, obj.pn_stimdata, obj.filelist.azidata{kk}, obj.pn_data, overwrite_flag);
                currdata = importdata(obj.filelist.azidata{kk});
                isVisuallyResponsive_azi(:, kk) = currdata.isVisuallyResponsive;
                RespMat_azi(:, :, :, :, kk) = currdata.RespVec;
            end
            obj.RFdata.RFmapping_results = structPacker([], isVisuallyResponsive_alt, 'isVisuallyResponsive_alt', isVisuallyResponsive_azi, 'isVisuallyResponsive_azi');
        end
        
        function results_struct = processorB(obj)
            alt_fit = cell(obj.num_sessions, obj.num_cells);
            azi_fit = cell(obj.num_sessions, obj.num_cells);
            alt_pref = zeros(obj.num_sessions, obj.num_cells);
            azi_pref = zeros(obj.num_sessions, obj.num_cells);
            alt_p = zeros(obj.num_sessions, obj.num_cells);
            azi_p = zeros(obj.num_sessions, obj.num_cells);
            isSpatiallyTuned = zeros(obj.num_sessions, obj.num_cells);
            isSpatiallyTunedAlt = zeros(obj.num_sessions, obj.num_cells);
            isSpatiallyTunedAzi = zeros(obj.num_sessions, obj.num_cells);
            roi_centroids = zeros(obj.num_cells, 2, obj.num_sessions);
            alt_corr = zeros(1, obj.num_sessions);
            azi_corr = zeros(1, obj.num_sessions);
            
            for kk = 1:obj.num_sessions
                fprintf('Analyzing session %d out of %d\n', kk, obj.num_sessions);               
                RFmapping_curr = ReceptiveField_analysis(obj.filelist.altdata{kk}, obj.pn_data, obj.filelist.azidata{kk}, obj.pn_data, 'No');      %Run sequentially through the multiselected data files using corresponding alt and azi data
                alt_fit(kk, :)           = RFmapping_curr.alt_fit;
                azi_fit(kk, :)           = RFmapping_curr.azi_fit;
                alt_pref(kk, :)          = RFmapping_curr.alt_pref;
                azi_pref(kk, :)          = RFmapping_curr.azi_pref;
                alt_p(kk, :)             = RFmapping_curr.alt_p;
                azi_p(kk, :)             = RFmapping_curr.azi_p;
                isSpatiallyTuned(kk, :)  = RFmapping_curr.isSpatiallyTuned;
                isSpatiallyTunedAlt(kk, :)  = RFmapping_curr.isSpatiallyTunedAlt;
                isSpatiallyTunedAzi(kk, :)  = RFmapping_curr.isSpatiallyTunedAzi;
                roi_centroids(:, :, kk)  = RFmapping_curr.roi_centroids;
                alt_corr(kk)             = RFmapping_curr.alt_corr;
                azi_corr(kk)             = RFmapping_curr.azi_corr;
            end
     
            
            results_struct = structPacker(obj.RFdata.RFmapping_results, alt_fit, 'alt_fit', azi_fit, 'azi_fit', alt_pref, 'alt_pref',...
                    azi_pref, 'azi_pref', alt_p, 'alt_p', azi_p, 'azi_p', isSpatiallyTuned, 'isSpatiallyTuned', isSpatiallyTunedAlt, 'isSpatiallyTunedAlt', isSpatiallyTunedAzi, 'isSpatiallyTunedAzi',...
                    roi_centroids, 'roi_centroids', alt_corr, 'alt_corr', azi_corr, 'azi_corr');
        end
        
        %Getter
        function RFdata = getRFdata(obj)
            RFdata = obj.RFdata;
        end
    end
end
        
        