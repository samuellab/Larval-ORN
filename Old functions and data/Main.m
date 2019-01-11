%Fig. 3, Sup Fig. 5 
PCAonData();            %PCA on dose-response data 

% Fig. 2
PlotHeatMap();          % visualize the raw data using heatmap

% Fig. 4C
Plot_AvsC();            % plot the activity vs concentration

% Fig. 4A 4B, Sub Fig. 5A, 6A
DoseResponseAnalysis(); % fit the does response curves

% Fig. 4D, Sub Fig 6B, 6C
AnalyzeEC50Matrix();    % analyze the EC50 matrix

% Fig. 6, Sup Fig.9
SimuPN(); % calculate uPN population activities using divisive normalizaiton model

% Fig. 6, Sup Fig.9
SimmPN(); %simulat the mPN circuit, using LIF model, given wiring diagram and ORN input