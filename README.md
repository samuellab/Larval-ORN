# Larval-ORN

Data and analysis of fruit fly larvae Olfactory Receptor Neurons (ORNs) responses to a broad panel of odors. This is the code and data from the paper, "Structured odorant response patterns across a complete olfactory receptor neuron population".



Quick Links
--------------------
* Raw dose-response data (Fig. 2A): Larval-ORN/Figure2/results/ORN_Dose_Response_Fig2A.xlsx
* Sensitivity matrix (Fig. 3B): Larval-ORN/Figure3/results/log_10_EC50.csv

Instructions
--------------------------------
All code was written using Matlab 2018b.
Code and data are organized by figure in order to ease finding and reproducing any figure from the paper.
Each "figure" folder includes a "data", "functions", and "results" subfolder, as well as matlab script(s) to generate panels within each figure. Below is a summary of the contents that can be found in each "figure" folder:

Figure 1 folder:
* AutoCAD design file of the microfluidic chips
* Parts list of materials used to build the odor delivery setup

Figure 2 folder:
* 21 larval ORNs' activity data in response to 34 odorants at five orders of magnitude in concentration
* Principal Component Analysis of ORN dose-response data
* PCA analysis of ORN dynamic data
* ORN calcium imaging response traces
* Analysis of our selected 35 odors sampled from odor space explored by the fly community

Figure 3 folder:
* Fit of the dose-response data using Hill equation
* Maximum Likelihood Estimation of ORN-odor sensitivies
* Fit of the sensitivity distribution
* Power-law relationship between average ORN activity and odor concentration
* Fit comparison to electrophysiological data

Figure 4 folder:
* Principal Component Analysis of the sensitivity matrix
* Comparision with odor molecular structure properties

Figure 5 folder:
* Reverse-correlation analysis on dynamic ORN activity responding to white-noise odor stimuli
* Comparision of linear filters acorss animals


How to cite:
----------------
Please cite the following paper if you use the data or code:

> Si G., Kanwal J.K., Hu Y., Tabone C., Baron J., Berck M., Vignoud G., Samuel A.D.T.(2019). Structured Odorant Response Patterns across a Complete Olfactory Receptor Neuron Population. Neuron, 101, 1-13.
