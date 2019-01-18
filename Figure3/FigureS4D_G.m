% each script can be run independently, all will save data to ./results

%% Perform MLE fitting to find distributions for A, x0, n

MLEScript; %very slow

%% Perform Fitting to find parameters for each trial

invMLEFit;


%% Plot distributions in Figure S4D-G

plotInvDistributions;
