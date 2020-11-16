#!/bin/bash

jobID="470"

#root -l -b -q plotAll.cxx\(\"${jobID}\"\)
#root -l -b -q m2BypTBins.cxx\(\"${jobID}\"\)
#root -l -b -q yVsEtaPlots.cxx\(\"${jobID}\"\)
root -l -b -q acceptanceCuts.cxx\(\"${jobID}\"\)

# Integrated yields
root -l -b -q intYield.cxx\(\"${jobID}\",\"E\"\)
root -l -b -q showAllFits.cxx\(\"${jobID}\",\"E\"\)
#root -l -b -q overlay.cxx\(\"${jobID}\",\"E\"\)

# EP Resolution and Flow Calculations
#root -l -b -q resolutions.cxx\(\"${jobID}\"\)
#root -l -b -q coefficients.cxx\(\"${jobID}\"\)
