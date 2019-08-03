# Larval-ORN

Data and analysis of fruit fly larvae olfactory receptor neuron (ORN) population responses to a broad panel of odors. This folder contains code and data from the paper, "Structured odorant response patterns across a complete olfactory receptor neuron population."



Quick Links
--------------------
* Raw dose-response data (Fig. 2A): Larval-ORN/Figure2/results/ORN_Dose_Response_Fig2A.xlsx
* Sensitivity matrix (Fig. 3B): Larval-ORN/Figure3/results/log_10_EC50.csv

Instructions
--------------------------------
All code was written using Matlab 2018b.
Code and data are organized by figure to ease locating and reproducing any of the figures from the paper.
Each "figure" folder contains a "data", "functions", and "results" subfolder, as well as Matlab script(s) to generate the subfigures. Below is a summary of the contents that can be found in each "figure" folder:

Figure 1 folder:
* AutoCAD design file of the microfluidic chips
* Parts list of materials used to build the odor delivery setup

Figure 2 folder:
* Activity data of the 21 larval ORNs in response to 34 odorants at five orders of magnitude in concentration
* Principal Component Analysis of ORN dose-response data
* PCA analysis of ORN response over time data
* ORN calcium imaging response traces
* Sampling analysis of our selected 35 odors in comparison to those previously explored by the fly community

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
* Comparision of linear filters across animals


How to cite:
----------------
Please cite the following in press paper if you use the data or code:

> Si G., Kanwal J.K., Hu Y., Tabone C., Baron J., Berck M., Vignoud G., Samuel A.D.T. Structured Odorant Response Patterns across a Complete Olfactory Receptor Neuron Population. Neuron 101.5 (2019): 950-962. https://doi.org/10.1016/j.neuron.2018.12.030
