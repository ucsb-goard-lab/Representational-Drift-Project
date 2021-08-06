classdef ProcessorController < handle
    
    properties
        stimtype string
        processor           %subclass processor
    end
    
    methods
        
        function obj = ProcessorController(varargin)
            if ~isempty(varargin)
                obj.stimtype = varargin(1);
            end
            switch obj.stimtype
                case 'PDG'
                    obj.processor = PDG_Processor();
                case 'NatMov'
                    obj.processor = NatMov_Processor();
                case 'RF'
                    obj.processor = RF_Processor();
                case 'MultiMov'
                    obj.processor = MultiMov_Processor();
                case 'MultiMovtoNatMov'
                    obj.processor = MultiMovtoNatMov_Processor();
                case 'NatScene'
                    obj.processor = NatScene_Processor();
                case 'NatMovdisc'
                    obj.processor = PDG_Processor();        % same exact structure as PDG
                case 'PDGcont'  
                    obj.processor = PDGcont_Processor();
            end
        end
        
        function runProcessor(obj)
            obj.processor.run();
        end
    end
end
                    