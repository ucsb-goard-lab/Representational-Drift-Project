classdef Eye_Processor < General_Processor

	properties
		processed_pupil       % processed eyetracking data struct
		eyeinfo               % raw imported eyetracking data
		eyemovie_stimdata     % eyetracking video descriptor statistics
		STIM_LIST = {'PDG', 'NatMov', 'MultiMov'};        % possible accompanying stimuli
	end

	properties (Access = protected)
		Resp                  % reference RespData imported for downsampling of eyetracking data
	end

	methods
		function obj = Eye_Processor()
    	end

    	function run(obj, tracker_flag)
    		obj.runA(tracker_flag);
    		obj.runB();
    	end

    	function runA(obj, tracker_flag)	
        % tracker_flag: 'Yes' or 'No', whether or not to run submethod A (video preprocessor)	
    		fprintf('Select movie stimdata files.\n');
			obj.importMoviedata();

			if strcmp(tracker_flag, 'Yes')
				obj.runReader();
				obj.runCleaner();
				obj.runCropper();
				obj.runPupilDetector('No');
			end
    	end

    	function runB(obj)
        % secondary processor, packages eyetracking data for saving
    		[Resp_referencefile, pn] = uigetfile('.mat');
    		cd(pn)
    		Resp_reference = importdata(Resp_referencefile);
            
            if obj.num_sessions == size(Resp_reference.PDG.RespMat_Full, 4)
                for ii = 1:obj.num_sessions
                    movieframes = length(obj.eyeinfo(ii).pupil);

                    for ff = 1:movieframes
                        raw_area(ff) = obj.eyeinfo(ii).pupil(ff).Area;
                        raw_centroid(ff, :) = obj.eyeinfo(ii).pupil(ff).Centroid;
                        raw_eccentricity(ff) = obj.eyeinfo(ii).pupil(ff).Eccentricity;
                        raw_orientation(ff) = obj.eyeinfo(ii).pupil(ff).Orientation;
                        raw_majoraxis(ff) = obj.eyeinfo(ii).pupil(ff).MajorAxisLength;
                        raw_minoraxis(ff) = obj.eyeinfo(ii).pupil(ff).MinorAxisLength;
                    end
                    area(ii).raw = raw_area;
                    centroid(ii).raw = raw_centroid;
                    eccentricity(ii).raw = raw_eccentricity;
                    orientation(ii).raw = raw_orientation;
                    majoraxis(ii).raw = raw_majoraxis;
                    minoraxis(ii).raw = raw_minoraxis;
                    currmovie_stimdata = obj.eyemovie_stimdata(ii);
                    stim_type = obj.eyemovie_stimdata(ii).stim_type;

                    obj.Resp = Resp_reference.(stim_type).RespMat_Full(:, :, :, ii);

                    newsize = size(obj.Resp, 1)*size(obj.Resp, 2);
                    area(ii).resampled = obj.downsampleEye(raw_area, newsize);
                    centroid(ii).resampled(:, 1) = obj.downsampleEye(raw_centroid(:, 1)', newsize);
                    centroid(ii).resampled(:, 2) = obj.downsampleEye(raw_centroid(:, 2)', newsize);
                    eccentricity(ii).resampled = obj.downsampleEye(raw_eccentricity, newsize);
                    orientation(ii).resampled = obj.downsampleEye(raw_orientation, newsize);
                    majoraxis(ii).resampled = obj.downsampleEye(raw_majoraxis, newsize);
                    minoraxis(ii).resampled = obj.downsampleEye(raw_minoraxis, newsize);

                    area(ii).sorted = obj.sortFrames(area(ii).resampled);
                    centroid(ii).sorted(:, :, 1) = obj.sortFrames(centroid(ii).resampled(:, 1));
                    centroid(ii).sorted(:, :, 2) = obj.sortFrames(centroid(ii).resampled(:, 2));
                    eccentricity(ii).sorted = obj.sortFrames(eccentricity(ii).resampled);
                    orientation(ii).sorted = obj.sortFrames(orientation(ii).resampled);
                    majoraxis(ii).sorted = obj.sortFrames(majoraxis(ii).resampled);
                    minoraxis(ii).sorted = obj.sortFrames(minoraxis(ii).resampled);

                    mean_frame(:, :, ii) = mean(obj.eyeinfo(ii).cropped_movie, 3);
                end
                processed_pupil.area = area;
                processed_pupil.centroid = centroid;
                processed_pupil.eccentricity = eccentricity;
                processed_pupil.orientation = orientation;
                processed_pupil.majoraxis = majoraxis;
                processed_pupil.minoraxis = minoraxis;
                processed_pupil.stim_type = stim_type;
                processed_pupil.movie_size = size(obj.eyeinfo(1).cropped_movie);
                processed_pupil.mean_frame = mean_frame;


                obj.processed_pupil = processed_pupil;
            elseif obj.num_sessions > size(Resp_reference.PDG.RespMat_Full, 4)
                for ii = 1:obj.num_sessions
                    movieframes = length(obj.eyeinfo(ii).pupil);
                    
                    raw_area = zeros(1, movieframes);
                    raw_centroid = zeros(movieframes, 2);
                    raw_eccentricity = zeros(1, movieframes);
                    raw_orientation = zeros(1, movieframes);
                    raw_majoraxis = zeros(1, movieframes);
                    raw_minoraxis = zeros(1, movieframes);
                    for ff = 1:movieframes
                        raw_area(ff) = obj.eyeinfo(ii).pupil(ff).Area;
                        raw_centroid(ff, :) = obj.eyeinfo(ii).pupil(ff).Centroid;
                        raw_eccentricity(ff) = obj.eyeinfo(ii).pupil(ff).Eccentricity;
                        raw_orientation(ff) = obj.eyeinfo(ii).pupil(ff).Orientation;
                        raw_majoraxis(ff) = obj.eyeinfo(ii).pupil(ff).MajorAxisLength;
                        raw_minoraxis(ff) = obj.eyeinfo(ii).pupil(ff).MinorAxisLength;
                    end
                    stim_type = obj.eyemovie_stimdata(ii).stim_type;
                    real_sessionidx = ceil(ii/2);
                    area(real_sessionidx).(stim_type).raw = raw_area;
                    centroid(real_sessionidx).(stim_type).raw = raw_centroid;
                    eccentricity(real_sessionidx).(stim_type).raw = raw_eccentricity;
                    orientation(real_sessionidx).(stim_type).raw = raw_orientation;
                    majoraxis(real_sessionidx).(stim_type).raw = raw_majoraxis;
                    minoraxis(real_sessionidx).(stim_type).raw = raw_minoraxis;
                    currmovie_stimdata = obj.eyemovie_stimdata(ii);

                    obj.Resp = Resp_reference.(stim_type).RespMat_Full(:, :, :, real_sessionidx);

                    newsize = size(obj.Resp, 1)*size(obj.Resp, 2);
                    area(real_sessionidx).(stim_type).resampled = obj.downsampleEye(raw_area, newsize);
                    centroid(real_sessionidx).(stim_type).resampled(:, 1) = obj.downsampleEye(raw_centroid(:, 1)', newsize);
                    centroid(real_sessionidx).(stim_type).resampled(:, 2) = obj.downsampleEye(raw_centroid(:, 2)', newsize);
                    eccentricity(real_sessionidx).(stim_type).resampled = obj.downsampleEye(raw_eccentricity, newsize);
                    orientation(real_sessionidx).(stim_type).resampled = obj.downsampleEye(raw_orientation, newsize);
                    majoraxis(real_sessionidx).(stim_type).resampled = obj.downsampleEye(raw_majoraxis, newsize);
                    minoraxis(real_sessionidx).(stim_type).resampled = obj.downsampleEye(raw_minoraxis, newsize);

                    area(real_sessionidx).(stim_type).sorted = obj.sortFrames(area(real_sessionidx).(stim_type).resampled);
                    centroid(real_sessionidx).(stim_type).sorted(:, :, 1) = obj.sortFrames(centroid(real_sessionidx).(stim_type).resampled(:, 1));
                    centroid(real_sessionidx).(stim_type).sorted(:, :, 2) = obj.sortFrames(centroid(real_sessionidx).(stim_type).resampled(:, 2));
                    eccentricity(real_sessionidx).(stim_type).sorted = obj.sortFrames(eccentricity(real_sessionidx).(stim_type).resampled);
                    orientation(real_sessionidx).(stim_type).sorted = obj.sortFrames(orientation(real_sessionidx).(stim_type).resampled);
                    majoraxis(real_sessionidx).(stim_type).sorted = obj.sortFrames(majoraxis(real_sessionidx).(stim_type).resampled);
                    minoraxis(real_sessionidx).(stim_type).sorted = obj.sortFrames(minoraxis(real_sessionidx).(stim_type).resampled);

                    mean_frame.(stim_type)(:, :, real_sessionidx) = mean(obj.eyeinfo(ii).cropped_movie, 3);
                end
                processed_pupil.area = area;
                processed_pupil.centroid = centroid;
                processed_pupil.eccentricity = eccentricity;
                processed_pupil.orientation = orientation;
                processed_pupil.majoraxis = majoraxis;
                processed_pupil.minoraxis = minoraxis;
                processed_pupil.stim_type = stim_type;
                processed_pupil.movie_size = size(obj.eyeinfo(1).cropped_movie);
                processed_pupil.mean_frame = mean_frame;


                obj.processed_pupil = processed_pupil;
            end
    	end

    	function runReader(obj)
        % imports eye tracking videos using readEyeTrackingVideo function
    		for ii = 1:obj.num_sessions
    			if ii == 1
    				obj.eyeinfo = EyeTracker();
    			else
    				obj.eyeinfo(ii) = EyeTracker();
    			end
	    		obj.eyeinfo(ii).readEyeTrackingVideo();
	    	end
    	end

    	function runCleaner(obj)
        % cleans up dropped frames 
    		for ii = 1:obj.num_sessions
	    		obj.eyeinfo(ii).cleanVideo('interpolate');
	    	end
	    end

	    function runCropper(obj)
        % crops all movies to common coordinate plane by resizing according to user-defined anchor points
	    	for ii = 1:obj.num_sessions
				obj.eyeinfo(ii).cropMovie('Points');
			end
			%%need to resize movies to reference cropped movie size
			if obj.num_sessions > 1		%skip if you only have 1 movie
				targetsize = size(obj.eyeinfo(1).cropped_movie);
				target_y = targetsize(1);
				target_x = targetsize(2);
				for ii = 2:obj.num_sessions
					curr_movie = obj.eyeinfo(ii).cropped_movie;
					new_movie = zeros(target_y, target_x, size(curr_movie, 3));
					for ff = 1:size(curr_movie, 3)
						new_movie(:, :, ff) = imresize(curr_movie(:, :, ff), [target_y target_x]);
					end
					% obj.eyeinfo(ii).orig_cropped_movie = curr_movie;		%preserve original cropped movie
					obj.eyeinfo(ii).cropped_movie = new_movie;		%replace cropped movie with newly resized movie
				end
			end
		end

		function runPupilDetector(obj, check)
        % tracks aspects of pupil across recording
			for ii = 1:obj.num_sessions
				obj.eyeinfo(ii).detectPupil();
				if strcmp(check, 'Yes')
					obj.eyeinfo(ii).checkPerformance(1:2000);
				end
			end
		end

    	function saveData(obj)
    		processed_pupil = obj.processed_pupil;
    		save pupilinfo processed_pupil

    		for ii = 1:obj.num_sessions
    			pupilraw{ii} = obj.eyeinfo(ii).pupil;
    		end
    		save pupilraw pupilraw
    	end

    	function importSave(obj)
    		[filename, pn] = uigetfile('.mat');
    		cd(pn);
    		pupilraw = importdata(filename);
    		obj.num_sessions = length(pupilraw);
    		for ii = 1:length(pupilraw)
    			if ii == 1
    				obj.eyeinfo = EyeTracker();
    			else
    				obj.eyeinfo(ii) = EyeTracker();
    			end
    			obj.eyeinfo(ii).pupil = pupilraw{ii};
    		end
    	end

    end

    methods (Access = private)

		function sorted_frames = sortFrames(obj, inputvec)
        % data sorter
			num_reps = size(obj.Resp, 1);
			frames_per_rep = size(obj.Resp, 2);
			sorted_frames = zeros(num_reps, frames_per_rep);
			for rep = 1:num_reps
				curr_frame = (rep-1)*frames_per_rep;
				sorted_frames(rep, :) = inputvec(curr_frame+1:curr_frame+frames_per_rep);
			end
		end

    	function importMoviedata(obj)
        % importing descriptive datafiles corresponding to eyetracking movies
    		[filenames, pn] = uigetfile('.mat', 'MultiSelect', 'on');
    		cd(pn);
    		obj.num_sessions = length(filenames);
    		for ii = 1:length(filenames)
    			temp = importdata(filenames{ii});
				for jj = 1:length(obj.STIM_LIST)
	    			curr_search = obj.STIM_LIST{jj};
	    			found = strfind(filenames{ii}, curr_search);
	    			if found
	    				temp.stim_type = curr_search;
	    				break
	    			end
	    		end
	    		eyemovie_stimdata(ii) = temp;
	    	end
	    	obj.eyemovie_stimdata = eyemovie_stimdata; 	%to avoid assignment issues
    	end

    	function out = downsampleEye(obj, inputvec, targetsize)
        % downsampling raw eyetracking metrics to match framerate of 2P data
    		rescale_factor = size(inputvec, 2)/targetsize;
    		out = zeros(1, targetsize);
    		index_tracker = 0;		%initialize tracker
    		for ff = 1:targetsize
    			nearest_idx = ceil(ff*rescale_factor);
    			num_frames_toavg = nearest_idx - index_tracker;
    			if index_tracker ~= 0
    				try
    					out(ff) = mean(inputvec(index_tracker:index_tracker+num_frames_toavg));
    				catch
    					out(ff) = mean(inputvec(index_tracker:index_tracker:end));
    				end
    			else
    				out(ff) = mean(inputvec(1:num_frames_toavg));
    			end
    			index_tracker = nearest_idx;
    		end
    	end
    end
end
