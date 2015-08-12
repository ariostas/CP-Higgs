# Files for Combine

## Use

To use this you must have combine installed. First modify add a file Histograms.root with the histograms, and modify cphiggs.txt appropriately. Then simply run the following commands.
```
python ../../CMSSW_7_1_5/src/HiggsAnalysis/CombinedLimit/scripts/text2workspace.py cphiggs.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:CPHiggsModel -o cphiggs.root
combine -M MaxLikelihoodFit cphiggs.root --plots
combine -M Asymptotic cphiggs.root
```
