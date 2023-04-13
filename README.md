# FlowAnalysis

# This program is deprecated and was replaced with TreeAnalsis

ROOT macro for calculating coefficients of anisotropic flow (of order_n) using the Event Plane Method (of order_m).

* This requires StRoot/StEpdUtil and StRoot/StPicoEvent libraries.

* The settings for the macro (FlowAnalyzer) are controlled by the config_\*.txt files. These are read and accessed via the ConfigReader.

* Q vector recentering, Fourier shifting, and event plane resolution corrections are implemented, so the analyzer has to be run 4 times and supplied the correctionInfo_INPUT.root (which it produces and updates automatically) and resolutionInfo_INPUT.root (which is produced in the postAnalysis) files during these iterations. The corrections info is supplied in the 2nd, 3rd, and 4th iterations, while the resolution info is only in the 4th.

* FlowAnalyzerABCD.cxx and FlowAnalyzerEF.cxx are deprecated.
