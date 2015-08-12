# CP phase measurement in H->TauTau decays

## Description

Here you can find the code used for the selection of H->TauTau decays in both Bacon and Delphes. It also includes the code for using Michael Williams' method and the CombinedLimit analysis.

## Use

To use this you simply need to go into the appropriate directory for Bacon or Delphes, and then run it with
```
root -b -q cphiggs.C+
```
This runs only the selection, however if you also want to run Mike's method you must use
```
root -b -q cphiggs.C+\(\"true\"\)
```
Note that this runs over small ntuples that were produces with the selection code found in the selection directory.

Running this produces a root file with the histograms that can later be used for Combine.