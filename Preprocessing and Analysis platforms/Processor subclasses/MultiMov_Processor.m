classdef MultiMov_Processor < ResponseProcessor & StabilityProcessor
    
    properties (Access = protected)
    end
    
    methods
        function obj = MultiMov_Processor()
        end
        
        function run(obj)
            obj.getFiles('MultiMov');
            
            obj.Stimdata = obj.getStimData('MultiMov3', obj.ref_session);
            
            RespMat_Full = obj.sortData;
            RespMat_onTime = obj.extractOnTime(RespMat_Full);
            
            [movieset_full, movieset_ontime] = obj.unweave(RespMat_Full, RespMat_onTime);
            obj.RespData = structPacker([], movieset_full, 'RespMat_Full', movieset_ontime, 'RespMat_onTime');
        
            numMovs = obj.Stimdata(1).interleaved;
            CCs = cell(1, numMovs);
            for ii = 1:numMovs
                CCs{ii} = obj.computeCCs(movieset_ontime{ii});
            end
            
            RDI = cell(1, numMovs);
            CC_ws = cell(1, numMovs);
            CC_bs = cell(1, numMovs);
            for jj = 1:numMovs
                [RDI{jj}, CC_ws{jj}, CC_bs{jj}, RDI_control{jj}] = obj.computeStability(movieset_ontime{jj});
            end
            
            obj.StabilityData = structPacker([], CCs, 'CCs', RDI, 'RDI', ...
                                            CC_ws, 'CC_ws', CC_bs, 'CC_bs', RDI_control, 'RDI_control');  
        end
        
        function [fullset, ontimeset] = unweave(obj, fullMat, ontimeMat)
            fprintf('Unweaving movies...\n')
            numMovs = obj.Stimdata(1).interleaved;
            fullset = cell(1, numMovs);
            ontimeset = cell(1, numMovs);
            for uu = 1:numMovs
                fullset{uu} = fullMat((uu-1)*obj.Stimdata(uu).movieReps+1:(uu-1)*obj.Stimdata(uu).movieReps+obj.Stimdata(uu).movieReps, :, :, :);
                ontimeset{uu} = ontimeMat((uu-1)*obj.Stimdata(uu).movieReps+1:(uu-1)*obj.Stimdata(uu).movieReps+obj.Stimdata(uu).movieReps, :, :, :);
            end
        end
    end
end