classdef EyeAnalyzer < StabilityAnalyzer
    
    properties
    	pupilinfo
    	pupilraw
    end
    
    methods

    	function obj = EyeAnalyzer()
    		%empty (use importData)
    	end

    	function density = drawCentroidMap(obj)
    		ysize = obj.pupilinfo.movie_size(1);
    		xsize = obj.pupilinfo.movie_size(2);
			ypts = 1:ysize;
			xpts = 1:xsize;
    		num_plots = ceil(obj.num_sessions/3);
    		%plot bounds
    		top = 18;
    		bottom = 82;
    		left = 40;
    		right = 120;

    		dmaps = figure;
    		frames = figure;

    		for ii = 1:obj.num_sessions
    			currpoints = obj.pupilinfo.centroid(ii).resampled;

    			density(:, :, ii) = histcounts2(currpoints(:, 2), currpoints(:, 1), ypts, xpts);
    			smoothmap = imgaussfilt(density(:, :, ii), 2);

				% outline = obj.getEyeOutline(ii);
				% smoothmap(outline) = max(max(smoothmap));

    			figure(dmaps)
				dm(ii) = subplot(num_plots, 3, ii);
				imagesc(xpts, ypts, smoothmap);
				% set(gca, 'XLim', xpts([left right]), 'YLim', ypts([top bottom]));
				set(gca, 'XLim', xpts([1 end]), 'YLim', ypts([1 end]));
				colormap parula
                pbaspect(dm(ii), [xsize ysize 1])
				hold on
				obj.drawAveragePupilBoundary(ii);

				figure(frames)
				fr(ii) = subplot(num_plots, 3, ii);
				% imagesc(obj.pupilinfo.mean_frame(top:bottom, left:right, ii));
				imagesc(obj.pupilinfo.mean_frame(:, :, ii));	%should normalize all of these		
				colormap gray
                pbaspect(fr(ii), [xsize ysize 1])
    		end

    		avgmap = mean(density, 3);
    		smoothmap = imgaussfilt(avgmap, 2);
    		figure
    		imagesc(smoothmap);
			% set(gca, 'XLim', xpts([left right]), 'YLim', ypts([top bottom]));
    		hold on
    		obj.drawAveragePupilBoundary(1);
    		mean_centroids = obj.getMeanCentroids();
            scatter(mean_centroids(:, 1), mean_centroids(:, 2));
            pbaspect(gca, [xsize ysize 1])
    	end

    end

	methods (Access = private)

		function drawAveragePupilBoundary(obj, session)
			center = mean(obj.pupilinfo.centroid(session).resampled, 1);		%average centroid location 
			major = mean(obj.pupilinfo.majoraxis(session).resampled);
			minor = mean(obj.pupilinfo.minoraxis(session).resampled);
			ori = mean(obj.pupilinfo.orientation(session).resampled);

			theta = 0 : 0.01 : 2*pi;
            x = minor/2 * cos(theta);
            y = major/2 * sin(theta);
            
            % rotation?
            R = [cosd(ori), -sind(ori);... % create rotation matrix
            sind(ori), cosd(ori)];
            
            rCoords = R * [x; y]; % apply transform
            
            xr = rCoords(1, :)';
            yr = rCoords(2, :)';
            
            plot(yr + center(1), xr + center(2), 'LineWidth', 2, 'Color', [1, 0, 0]);
        end

        function outline = getEyeOutline(obj, session)
        	currimg = obj.pupilinfo.mean_frame(:, :, session);
        	dark = currimg < min(min(currimg)) + 3.5*std(std(currimg, [], 1), [], 2);
        	outline = bwmorph(dark, 'remove');
        end

        function centroids = getMeanCentroids(obj, input)
            if nargin < 2
                input = obj.pupilinfo.centroid;
            end
            for jj = 1:length(input)
                curr_mean = mean(input(jj).resampled, 1);
                if isempty(curr_mean)
                    centroids(jj, :) = NaN;
                else
                    centroids(jj, :) = curr_mean;
                end
            end
        end

        function stds = getStdCentroids(obj)
        	for jj = 1:obj.num_sessions
        		stds(jj, 1) = std(obj.pupilinfo.centroid(jj).resampled(:, 1), [], 1);
        		stds(jj, 2) = std(obj.pupilinfo.centroid(jj).resampled(:, 2), [], 1);
        	end
        end

        function size = getPupilSize(obj)
            for kk = 1:obj.num_sessions
                size{kk, :} = obj.pupilinfo.area(kk).raw;
            end
        end
    end
end