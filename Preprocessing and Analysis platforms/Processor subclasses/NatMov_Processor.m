classdef NatMov_Processor < ResponseProcessor & StabilityProcessor
    
    properties (Access = protected)
    end
    
    methods
        
        function obj = NatMov_Processor()
        end
        
        function run(obj)
            obj.getFiles('NatMov')
            
            obj.Stimdata = obj.getStimData('NatMov', obj.ref_session);

            RespMat_Full = obj.sortData;        %these internal methods will query for dataType, and will use DFF or reconstructed DFF accordingly
            RespMat_onTime = obj.extractOnTime(RespMat_Full);
            CCs = obj.computeCCs(RespMat_onTime);

            [RDI, CC_ws, CC_bs, RDI_control] = obj.computeStability(RespMat_onTime);

            obj.RespData = structPacker(obj.RespData, RespMat_Full, 'RespMat_Full', RespMat_onTime, 'RespMat_onTime');
            obj.StabilityData = structPacker(obj.StabilityData, CCs, 'CCs', RDI, 'RDI', CC_ws, 'CC_ws', CC_bs, 'CC_bs', RDI_control, 'RDI_control');
        end
        
    end
end
        