classdef RFAnalyzer < Analyzer
    
    properties (Access = protected)
        RF_matrix           %Receptive field matricies
    end
    
    methods
        
        function obj = RFAnalyzer()
        end

        function [A, B, C] = RFoverTime(obj, constraint)
            %constraint indicates minimum number of sessions for which a cell must be spatially tuned to be included
            %in the analysis. In any case, cells must also be spatially tuned for the reference session (day 0).
            alt_size = 100;         % degrees in alt direction
            azi_size = 130;         % degrees in azi direction
            alt_deg = 3.3;          % degrees per alt location
            azi_deg = 3.25;          % degrees per azi location

            alt_pref = obj.RFdata.RFmapping_results.alt_pref*alt_deg;
            azi_pref = obj.RFdata.RFmapping_results.azi_pref*azi_deg;
            isSpatiallyTunedAlt = logical(obj.RFdata.RFmapping_results.isSpatiallyTunedAlt);
            isSpatiallyTunedAzi = logical(obj.RFdata.RFmapping_results.isSpatiallyTunedAzi);
            qual = obj.getUse_cells;
            sessions = obj.getUse_sessions;
            present = sum(sessions, 2) == obj.num_sessions;
            good_cells = qual & present';

            cellList_alt = find(isSpatiallyTunedAlt(1, :) & sum(isSpatiallyTunedAlt, 1) >= constraint & good_cells);
            cellList_azi = find(isSpatiallyTunedAzi(1, :) & sum(isSpatiallyTunedAzi, 1) >= constraint & good_cells);
            RDIalt = cell(1,length(cellList_alt));
            RDIazi = cell(1,length(cellList_azi));
            delta_alt = cell(1, length(cellList_alt));
            delta_azi = cell(1, length(cellList_azi));
            for jj = 1:length(cellList_alt)
                currCell = cellList_alt(jj);
                RDIalt{jj} = obj.StabilityData.NatMov.RDI(isSpatiallyTunedAlt(2:obj.num_sessions, currCell), currCell);
                delta_alt{jj} = abs(alt_pref(isSpatiallyTunedAlt(2:obj.num_sessions, currCell), currCell) - alt_pref(1, currCell));
            end
            for jj = 1:length(cellList_azi)
                currCell = cellList_azi(jj);
                RDIazi{jj} = obj.StabilityData.NatMov.RDI(isSpatiallyTunedAzi(2:obj.num_sessions, currCell), currCell);
                delta_azi{jj} = abs(azi_pref(isSpatiallyTunedAzi(2:obj.num_sessions, currCell), currCell) - azi_pref(1, currCell));
            end
            RDIalt_mean = cellfun(@mean, RDIalt);
            RDIazi_mean = cellfun(@mean, RDIazi);
            % RDIalt_SEM = cellfun(@std, RDI);
            %RDI_SEM = cellfun(@std, RDI)./sqrt(cellfun(@length, RDI));
            delta_alt_mean = cellfun(@mean, delta_alt);
            delta_alt_SEM = cellfun(@std, delta_alt)./sqrt(cellfun(@length, delta_alt));
            delta_azi_mean = cellfun(@mean, delta_azi);
            delta_azi_SEM = cellfun(@std, delta_azi)./sqrt(cellfun(@length, delta_azi));
                        
            bySession_delta_alt = zeros(1, obj.num_sessions-1);
            bySession_delta_azi = zeros(1, obj.num_sessions-1);
            bySession_delta_alt_SEM = zeros(1, obj.num_sessions-1);
            bySession_delta_azi_SEM = zeros(1, obj.num_sessions-1);
            for yy = 2:obj.num_sessions
                % tmp_alt = [];
                % tmp_azi = [];
                for ww = 1:length(cellList_alt)
                    currCell = cellList_alt(ww);
                    tmp_alt(ww) = alt_pref(yy, currCell) - alt_pref(1, currCell);
                end
                for ww = 1:length(cellList_azi)
                    currCell = cellList_azi(ww);
                    tmp_azi(ww) = azi_pref(yy, currCell) - azi_pref(1, currCell);
                end
                bySession_delta_alt(yy-1) = mean(tmp_alt);
                bySession_delta_alt_SEM(yy-1) = std(tmp_alt)/sqrt(length(tmp_alt));
                bySession_delta_azi(yy-1) = mean(tmp_azi);
                bySession_delta_azi_SEM(yy-1) = std(tmp_azi)/sqrt(length(tmp_azi));
                clear tmp_alt tmp_azi
            end

            % alt_avg = squeeze(mean(mean(obj.RFdata.RespData.altitude, 2), 4));         %average over frames and reps
            % azi_avg = squeeze(mean(mean(obj.RFdata.RespData.azimuth, 2), 4));         
            
            % for kk = 1:length(cellList)
            %     figure
            %     set(gcf, 'Position', [300 400 1000 400])
            %     currCell = cellList(kk);
            %     subplot(1, 2, 1)
            %     hold on
            %     for ss = 1:obj.num_sessions
            %         plot(alt_avg(currCell, :, ss))
            %         title(['pref alt ' num2str(alt_pref(ss, currCell)) ' | delta alt ' num2str(alt_pref(ss, currCell) - alt_pref(1, currCell)) ' | spat tuned ' num2str(isSpatiallyTuned(ss, currCell))])
            %         pause
            %     end
            %     hold off
            %     subplot(1, 2, 2)
            %     hold on
            %     for ss = 1:obj.num_sessions
            %         plot(azi_avg(currCell, :, ss))
            %         title(['pref azi ' num2str(azi_pref(ss, currCell)) ' | delta azi ' num2str(azi_pref(ss, currCell) - azi_pref(1, currCell)) ' | spat tuned ' num2str(isSpatiallyTuned(ss, currCell))])
            %         pause
            %     end
            %     hold off
            %     pause
            %     close
            % end


            % figure
            % scatter(delta_alt_mean, RDI_mean, 'filled');
            % hold on
            % errorbar(delta_alt_mean, RDI_mean, delta_alt_SEM, 'horizontal', 'LineStyle', 'none');
            % % errorbar(delta_alt_mean, RDI_mean, RDI_SEM, 'vertical', 'LineStyle', 'none');
            % % [m, gof] = fit(delta_alt_mean', RDI_mean', 'poly1');
            % xlim([-30 30])
            % % plot(m, delta_alt_mean, RDI_mean)
            % title(sprintf('n = %d', length(cellList)));
            % xlabel('Mean altitude change across sessions')
            % ylabel('Mean RDI across sessions')
            % axis square
            % hold off
            
            x_limit = 15;

            figure
            hold on
            for pp = 1:length(cellList_alt)
                delta_alt{pp}(delta_alt{pp} > x_limit) = x_limit;
                scatter(delta_alt{pp}, RDIalt{pp}, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
            end
            errorbar(delta_alt_mean, RDIalt_mean, delta_alt_SEM, 'horizontal', 'LineStyle', 'none');
            scatter(delta_alt_mean, RDIalt_mean, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
            [r, p] = corr(delta_alt_mean', RDIalt_mean');
            title(sprintf('n = %d, r = %0.2f, p = %0.2e', length(cellList_alt), r, p));
            xlabel('delta altitude (degrees)')
            ylabel('respective RDI')
            xlim([-1 x_limit])
            axis square
            box off
            hold off
            
            figure
            hold on
            for pp = 1:length(cellList_azi)
                scatter(delta_azi{pp}, RDIazi{pp}, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
            end
            errorbar(delta_azi_mean, RDIazi_mean, delta_azi_SEM, 'horizontal', 'LineStyle', 'none');
            scatter(delta_azi_mean, RDIazi_mean, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
            [r, p] = corr(delta_azi_mean', RDIazi_mean');
            title(sprintf('n = %d, r = %0.2f, p = %0.2e', length(cellList_azi), r, p));
            xlabel('delta azimuth (degrees)')
            ylabel('respective RDI')
            % xlim([-5 30])
            axis square
            box off
            hold off
            
            figure
            errorbar(2:obj.num_sessions, bySession_delta_alt, bySession_delta_alt_SEM)
            hold on
            errorbar(2:obj.num_sessions, bySession_delta_azi, bySession_delta_azi_SEM)
            legend alt azi
            xlim([1 obj.num_sessions+1])
            ylim([-30 30])
            axis square
            box off

            A.azi.curve = bySession_delta_azi;
            A.azi.errorbars = bySession_delta_azi_SEM;
            A.alt.curve = bySession_delta_alt;
            A.alt.errorbars = bySession_delta_alt_SEM;

            B.gray.x = delta_alt;  B.gray.y = RDIalt;
            B.blue.x = delta_alt_mean;  B.blue.y = RDIalt_mean;  B.blue.errorbars = delta_alt_SEM;

            C.gray.x = delta_azi;  C.gray.y = RDIazi;
            C.blue.x = delta_azi_mean;  C.blue.y = RDIazi_mean;  C.blue.errorbars = delta_azi_SEM;
        end
        
        function RF_matrix = extractRFs(obj)
            mean_altitude = squeeze(mean(mean(obj.RespData.RF.altitude, 2), 4));          %FOR V2: response matrices now exist in RespData, not RFmapping_results, and dimension order is different
            mean_azimuth = squeeze(mean(mean(obj.RespData.RF.azimuth, 2), 4));
            isVisuallyResponsive = (obj.RFdata.RFmapping_results.isVisuallyResponsive_alt | obj.RFdata.RFmapping_results.isVisuallyResponsive_azi);
            RF_matrix = zeros(obj.numCells, size(mean_azimuth, 2), size(mean_altitude, 2), obj.num_sessions);
            for rr = 1:size(mean_altitude, 3)
                for ii = 1:size(isVisuallyResponsive, 1)
                    hor = mean_altitude(ii,:, rr);
                    vert = mean_azimuth(ii,:, rr);
                    RF_matrix(ii,:,:, rr) = bsxfun(@times,hor,vert');
                end    
            end
        end
        
        %Getter
        function RF_matrix = getRFmatrix(obj)
            RF_matrix = obj.RF_matrix;
        end
    end
end