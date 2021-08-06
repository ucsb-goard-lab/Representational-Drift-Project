
The 'Preprocessing and Analysis platforms' folder contains software for turning .mat files from the intial pipeline ('Goard Lab 2P pre-processing pipeline') 
into further processed and more organized data structures to be used for analysis. These processed data structures are the format in which the data for this study is provided. 
The 'Stability_masterscript' script demonstrates this process. In addition to the classes in this folder that are responsible for generating the data in the format provided,
there are also classes used in analysis and the generation of figures by the scripts noted below. 

The 'Main figures' folder contains scripts for the analyses performed and figures generated for the main figures. 
The 'Other figures' folder contains scripts and functions for analysis done for supplementary figures and otherwise. 
(much of this code utilizes code from the folder mentioned above). 


The 'helper functions' folder contains various functions utilized by the code in the other folders.


The 'Goard Lab 2P pre-processing pipeline' folder contains software for processing 2-photon calcium imaging data in the form of multi-tif files. 
This pipeline (A > B > C) produces .mat files containing neural response DFF data in [neurons x frames] format. The versions used for this study are included here. 
An updated version of this code is available at: https://github.com/ucsb-goard-lab/Two-photon-calcium-post-processing
This folder also includes the software used for producing deconvolved spike data from the DF/F data (used for event detection in Figure 2).



Notes on the formatting of the figure generation code:

Unless explicitly formatted as a full function, the scripts for MS analyses and figure generation are meant to be run in a piecewise manner depending on the analysis
of interest, instead of all at once. Analyses are titled and separated into MATLAB-formatted code sections (denoted by %%). For guidance, refer to the comments. 
The analysis and figure generation scripts frequently use the provided Analysis platform, which utilizes Object Oriented Programming (OOP). Objects contain 
properties and methods, both of which are utilized using dot notation. For example:

all = StabilityAnalyzer; 	% A StabilityAnalyzer 'object' is created by generating an instance 
				% of the StabilityAnalyzer class.
all.importData;			% Associated with this class is a number of 'methods' that act as functions 
				% that can be called using dot notation. The importData method prompts you 
				% to select and import the datafiles necessary for running other methods for 
				% this class. When importing data in this way for any of the classes in the
				% package, it is necessary to simultaneously select the following datafiles:
				% 1) RespData.mat 2) RoiINFO.mat 3) StabilityData.mat
				% If pooling data using the addData method (see below), all of these datafiles 
				% must be imported for each subsequent imaging field added.
				% The OrientationData.mat datafile should also be imported 
				% along with the above if you intend to perform any of the methods or analyses 
				% that require it. 
all.addData;			% The first imaging field loaded into the object must be performed with the 
				% importData method, but each subsequent imaging field to be pooled must be added 
				% using the addData method.
all.addData([1 1 0 1 1 1]);	% If a mouse being added has any recordings which occur on non-consecutive weeks, 
				% it must be added by using an input argument in the form of a vector denoting the 
				% weeks to which the data corresponds. The example here denotes data collected 
				% on weeks 1, 2, 4, 5, 6. 
all.setQualityThreshold(3); 	% This method sets the ROI quality threshold for inclusion in subsequent analyses.
all.cellSelection; 		% This method prompts the user to select the criteria by which ROIs are included 
				% in analyses. For instance, simultaneously select 'quality', 'presence', 
				% 'PDG_Responsive_thresh', and 'NatMov_Responsive_thresh' to include neurons 
				% above the set quality threshold and responsive to both the PDG and MOV (NatMov) stimuli.
					
				% Review each class's code for a full list of its methods
