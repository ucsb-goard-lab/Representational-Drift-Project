

JCS003 = RFAnalyzer();
JCS003.importData();
JCS003.setQualityThreshold(3);
JCS003.cellSelection;
[A, B, C] = JCS003.RFoverTime(7);           % input arg: minimum number of sessions for which a cell must be spatially tuned to be included