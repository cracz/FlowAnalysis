#!/bin/bash


jobID="26F"
order_n="3"

#root -l -b -q plotAll.cxx\(\"${jobID}\"\)
#root -l -b -q m2BypTBins.cxx\(\"${jobID}\"\)
#root -l -b -q yVsEtaPlots.cxx\(\"${jobID}\"\)
#root -l -b -q acceptanceCuts.cxx\(\"${jobID}\"\)
#root -l -b -q eventPlanes.cxx\(\"${jobID}\"\)
#root -l -b -q correlations.cxx\(\"${jobID}\",\"${order_n}\"\)

# Integrated yields
#root -l -b -q intYield.cxx\(\"${jobID}\",\"E\"\)
#root -l -b -q showAllFits.cxx\(\"${jobID}\",\"E\"\)
#root -l -b -q overlay.cxx\(\"${jobID}\",\"E\"\)

# EP Resolution and Flow Calculations
#root -l -b -q resolutions.cxx\(\"${jobID}\",\"${order_n}\"\)
root -l -b -q coefficients.cxx\(\"${jobID}\",\"${order_n}\"\)
#root -l -b -q vnVsY.cxx\(\"${jobID}\",\"${order_n}\"\)
#root -l -b -q vnVsPt.cxx\(\"${jobID}\",\"${order_n}\"\)
#root -l -b -q vnVsKt.cxx\(\"${jobID}\",\"${order_n}\"\)
#root -l -b -q finalWithSystematics.cxx\(\"${order_n}\"\)
#root -l -b -q meanpT.cxx\(\"${jobID}\",\"${order_n}\"\)
