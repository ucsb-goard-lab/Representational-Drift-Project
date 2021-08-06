classdef MultiMovtoNatMov_Processor < handle
    %special class that instantiates the MultiMov processor class, doesn't need to inherit from anything
    % processes Scrambled movie data as normal, but extracts normal MOV trials and saves the data on its own

    % RECOMMENDED TO USE THE OTHER SCRIPT FOR THIS, just takes the MOV data directly from saved scrambled movie datafiles


    properties (Access = protected)
        natmov
    end
    
    methods 
        function obj = MultiMovtoNatMov_Processor()
            obj.natmov = ProcessorController('natmov'); %create instance of natmov class
        end

        function run(obj)
            init_multimov_processor = MultiMov_Processor();
            init_multimov_processor.setCCmethod('Random');
            init_multimov_processor.setSubsampleflag('Yes');
            init_multimov_processor.setDataType('DFF');
            init_multimov_processor.run();

            multimov_RespData = init_multimov_processor.getRespData;
            multimov_StabilityData = init_multimov_processor.getStabilityData;

            obj.extractNatMovInfo(multimov_RespData, multimov_StabilityData); 
        end

        function extractNatMovInfo(obj, origRespData, origStabilityData)
            extractedRespData = structPacker([], origRespData.RespMat_Full{1}, 'RespMat_Full', origRespData.RespMat_onTime{1}, 'RespMat_onTime');
            extractedStabilityData = structPacker([], origStabilityData.MultiMov_CCs{1}, 'NatMov_CCs', origStabilityData.MultiMov_RDI{1}, 'NatMov_RDI',...
                origStabilityData.MultiMov.CC_ws{1}, 'CC_ws', origStabilityData.MultiMov.CC_bs{1}, 'CC_bs', origStabilityData.MultiMov.RDI_control{1}, 'RDI_control');
            obj.natmov.setRespData(extractedRespData);
            obj.natmov.setStabilityData(extractedStabilityData);
        end
        
        %Getter
        function out = getNatMov(obj)
            out = obj.natmov;
        end
    end
end
      
    