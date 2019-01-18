# Larval-ORN

Data and analysis of fruit fly larvae Olfactory Receptor Neruons (ORNs) activity responding to a broad panel of odors. It is the code and data part of the paper "Structured odorant response patterns across a complete olfactory receptor neuron population".



Quick Links
--------------------
Links to the following data sets:
1. Raw dose-response data (Fig. 2A): https://github.com/samuellab/Larval-ORN/tree/master/Figure2/Data
2. Averaged dose-response data of saturated curves (Fig. 3A): https://github.com/samuellab/Larval-ORN/tree/master/Figure3/Data
3. Sensitivity Matrix (Fig. 3B): https://github.com/samuellab/Larval-ORN/tree/master/Figure3/Data

Links to the following methods:
1. MLE (Fig. 3): https://github.com/samuellab/Larval-ORN/tree/master/Figure3/Script
2. Curve-fitting (Fig. 3): https://github.com/samuellab/Larval-ORN/tree/master/Figure3/Script
3. Linear-Nonlinear model (Fig. 5): https://github.com/samuellab/Larval-ORN/tree/master/Figure5/Script
(csv file)


Instructions
--------------------------------
All code was written using Matlab 2018b.
Codes and data are organized to reproduce each figures in the paper.
Each figure folder include script to generate subfigures and data, function and result subfolders.

Figure 1 folder:
* AutoCAD design file of the microfluidic chips.
* Parts list of materials used to build the odor delivery setup

Figure 2 folder:
* 21 larval ORNs' activity data in responds to 34 odors at five order of magnitude concentration levels.
* Principal Component Analysis of ORN dose-response data
* PCA analysis of ORN dynamic data.
* ORN response traces.
* The selected 35 odors span the odor space explored by the fly community.

Figure 3 folder:
* Fit the dose-response data using Hill equation
* Maximum Likelihood Estimation of the ORN-odor sensitivies.
* Fit the sensitivity distribution.
* Power-law relationship between average ORN activity and odor concentration.
* Comparision with electrophysiological data

Figure 4 folder:
* Principal Component Analysis of the sensitivity matrix.
* Comparision with odor molecular structure properties.

Figure 5:
* Reverse-correlation analysis ORN's dynamic activity responding to white-noise odor stimuli.
* Comparision of linear filters acorss animals.


How to cite:
----------------
Please cite the following paper if you use the data or code:

> Si G., Kanwal J.K., Hu Y., Tabone C., Baron J., Berck M., Vignoud G., Samuel A.D.T.(2019). Structured Odorant Response Patterns across a Complete Olfactory Receptor Neuron Population. Neuron, 101, 1-13.
