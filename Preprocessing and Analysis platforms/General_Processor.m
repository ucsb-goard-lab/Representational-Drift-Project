classdef General_Processor < handle
  
    properties
        num_sessions        % number of sessions for the field
        num_cells           % number of cells for this field
        filelist            % list of imported files
        reference_number    % reference session designation (defaults to 1st session)
        ref_session         % data from reference session (for getting stimulus data)
        Stimdata            % Stimulus information for sorting data
    end
    
    methods       
        function obj = General_Processor()
            obj.Stimdata = struct();
        end

        function out_stimdat = get(obj, stimtype, ref_session)
        % fetches stimulus data for sorting DFF data
            if nargin < 3
                fprintf('No reference session yet imported, select recording file for framerate reference.\n');
                [temp, pn] = uigetfile('.mat');
                cd(pn);
                ref_session = importdata(temp);
            end
            switch stimtype
                case 'PDG'            %set to standard PDG timing
                    out_stimdat.numReps = 8;
                    out_stimdat.presTime = 2;
                    out_stimdat.interoffTime = 4;
                    out_stimdat.setsperRepeat = 12;
                    out_stimdat.offTime = 0;
                    out_stimdat.onTime = (out_stimdat.presTime + out_stimdat.interoffTime)*out_stimdat.setsperRepeat;
                    out_stimdat.repeatTime = out_stimdat.onTime + out_stimdat.offTime;
                    out_stimdat.framerate = ref_session.frameRate;
                    out_stimdat.interleaved = 1;           %how many sets of this stimuli are in the recording (to be sorted)
                    out_stimdat.totalReps = out_stimdat.numReps*out_stimdat.interleaved;
                case 'NatMov'           %set to standard NatMovToE timing
                    out_stimdat.numReps = 30;
                    out_stimdat.presTime = 30;
                    out_stimdat.interoffTime = 0;
                    out_stimdat.setsperRepeat = 1;
                    out_stimdat.offTime = 5;
                    out_stimdat.onTime = (out_stimdat.presTime + out_stimdat.interoffTime)*out_stimdat.setsperRepeat;
                    out_stimdat.repeatTime = out_stimdat.onTime + out_stimdat.offTime;
                    out_stimdat.framerate = ref_session.frameRate;       
                    out_stimdat.interleaved = 1;
                    out_stimdat.totalReps = out_stimdat.numReps*out_stimdat.interleaved;
                case 'MultiMov3'                 %3 scrambled movies (add another option for different amounts, etc)
                    disp('Select stimulus data files.')
                    stimulusfilenames = uigetfile('.mat', 'MultiSelect', 'on');
                    out_stimdat = struct;
                    for ii = 1:length(stimulusfilenames)
                        stimulusdata_curr = importdata(stimulusfilenames{ii});

                        out_stimdat(ii).presTime = 30;
                        out_stimdat(ii).interoffTime = 0;
                        out_stimdat(ii).setsperRepeat = 1;
                        out_stimdat(ii).offTime = 5;
                        out_stimdat(ii).onTime = (out_stimdat(ii).presTime + out_stimdat(ii).interoffTime)*out_stimdat(ii).setsperRepeat;
                        out_stimdat(ii).repeatTime = out_stimdat(ii).onTime + out_stimdat(ii).offTime;
                        out_stimdat(ii).movieReps = stimulusdata_curr.repeats;
                        out_stimdat(ii).framerate = ref_session.frameRate;       
                        out_stimdat(ii).interleaved = stimulusdata_curr.numMovs;
                        out_stimdat(ii).totalReps = out_stimdat(ii).movieReps*out_stimdat(ii).interleaved;
                        out_stimdat(ii).trialOrder = stimulusdata_curr.trialOrder;              %trials go from least scrambled (1) to most scrambled (n)
                    end
                case 'RF'          %note, these names have different definitions for altitude and azimuth stim information, to be used in separate preprocessor function
                    out_stimdat.numReps = 10;
                    out_stimdat.altlocs = 30;
                    out_stimdat.azilocs = 40;
                    out_stimdat.onTime = 1;
                    out_stimdat.framerate = ref_session.frameRate;

                case 'PDGcont'
                    out_stimdat.numReps = 30;
                    out_stimdat.presTime = 30;
                    out_stimdat.interoffTime = 0;
                    out_stimdat.setsperRepeat = 1;
                    out_stimdat.offTime = 4;
                    out_stimdat.onTime = (out_stimdat.presTime + out_stimdat.interoffTime)*out_stimdat.setsperRepeat;
                    out_stimdat.repeatTime = out_stimdat.onTime + out_stimdat.offTime;
                    out_stimdat.framerate = ref_session.frameRate;       
                    out_stimdat.interleaved = 1;
                    out_stimdat.totalReps = out_stimdat.numReps*out_stimdat.interleaved;
               case 'NatMovdisc'        % same structure as PDG
                    out_stimdat.numReps = 8;
                    out_stimdat.presTime = 2;
                    out_stimdat.interoffTime = 4;
                    out_stimdat.setsperRepeat = 12;
                    out_stimdat.offTime = 0;
                    out_stimdat.onTime = (out_stimdat.presTime + out_stimdat.interoffTime)*out_stimdat.setsperRepeat;
                    out_stimdat.repeatTime = out_stimdat.onTime + out_stimdat.offTime;
                    out_stimdat.framerate = ref_session.frameRate;
                    out_stimdat.interleaved = 1;           %how many sets of this stimuli are in the recording (to be sorted)
                    out_stimdat.totalReps = out_stimdat.numReps*out_stimdat.interleaved;

            end
        end
        
        function getFiles(obj, stimtype)
            sprintf('Select all %s files', stimtype);
            [obj.filelist, pathname] = uigetfile('.mat', 'MultiSelect', 'on');
            cd(pathname);
            if iscell(obj.filelist)
                lengthList = length(obj.filelist);
            else
                lengthList = 1;
            end
            
            if isempty(obj.reference_number)
                obj.reference_number = 1;          %reference session defaults to 1
            end

            obj.num_sessions = lengthList;
            if iscell(obj.filelist)
                obj.ref_session = importdata(obj.filelist{obj.reference_number});
            else
                obj.ref_session = importdata(obj.filelist);
            end
            obj.num_cells = length(obj.ref_session.cellMasks);
        end     
        
    end
end
        
