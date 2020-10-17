#!/bin/bash

jobID="shaoweiEPD_etaSeparated"

#root -l -b -q subsEFpostAnalysis.cxx\(\"${jobID}\"\)
#root -l -b -q plotAll.cxx\(\"${jobID}\"\)
root -l -b -q m2BypTBins.cxx\(\"${jobID}\"\)
#root -l -b -q yVsEtaPlots.cxx\(\"${jobID}\"\)
#root -l -b -q resolutions.cxx\(\"${jobID}\"\)
#root -l -b -q acceptanceCuts.cxx\(\"${jobID}\"\)
