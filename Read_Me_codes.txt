You need the following toolboxes for preprocessing:
1. nirs-toolbox
2. NIRS-GUI

You will need the following toolboxes specific to functional community analysis:
1. Network Community Toolbox (Simply add the toolbox to the MATLAB path)
2. Brain Connectivity Toolbox (Simply add the toolbox to the MATLAB path)
3. GenLouvain (To install the toolbox follow the instructions here: https://github.com/GenLouvain/GenLouvain)
_____________________________________________
MATLAB codes:
1. preprocessing_steps.m : Includes preprocessing steps & the function to save results for next steps
2. flex_alleg_calculate.m : Calculates flexibility(switching rate) and Allegiance matrix (It also calculates the features for some null networks that were not used in this study.)
3. dynFC_analysis.m : Does the analysis included in the manuscript (it also includes some analyses were not used in the study.)