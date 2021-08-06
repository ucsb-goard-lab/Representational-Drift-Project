classdef PDGcont_Processor < ResponseProcessor & StabilityProcessor
    
    properties (Access = protected)
    end
    
    methods
        
        function obj = PDGcont_Processor()
        end
        
        function run(obj)
            obj.getFiles('PDGcont')
            
            obj.Stimdata = obj.getStimData('PDGcont', obj.ref_session);

            RespMat_Full = obj.sortData;        %these internal methods will query for dataType, and will use DFF or reconstructed DFF accordingly
            RespMat_onTime = obj.extractOnTime(RespMat_Full);
            CCs = obj.computeCCs(RespMat_onTime);

            [RDI, CC_ws, CC_bs, RDI_control] = obj.computeStability(RespMat_onTime);

            obj.RespData = structPacker(obj.RespData, RespMat_Full, 'RespMat_Full', RespMat_onTime, 'RespMat_onTime');
            obj.StabilityData = structPacker(obj.StabilityData, CCs, 'CCs', RDI, 'RDI', CC_ws, 'CC_ws', CC_bs, 'CC_bs', RDI_control, 'RDI_control');
        end
        
    end
end