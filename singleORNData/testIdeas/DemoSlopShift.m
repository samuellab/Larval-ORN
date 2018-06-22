%% input info
input.matFiles = '22c_data.mat';
input.varNames = 'pp';
input.ORNs = 'Or22ca';
input.odors = {'methyl salicylate', 'anisole'};
    
input.concList = {[10^-11; 3.16*10^-11; 10^-10; 3.16*10^-10; 10^-9; 3.16*10^-9; 10^-8; 3.16*10^-8], ...
    [3.16*10^-8; 10^-7; 3.16*10^-7; 10^-6; 3.16*10^-6; 10^-5; 3.16*10^-5; 10^-4];}; % concentation list

% other settings
cColor =[0 0.4470 0.7410; 0.85 0.325 0.0980];

%% load data and setup fitting
load(input.matFiles, input.varNames);   % load files

if  strcmp(input.varNames, 'pp')
    dataPool = pp';
elseif strcmp(input.varNames, 'sigMat')
    dataPool = sigMat;
else
    error('Not listed input variable.');
end

dataOdor1 = dataPool(:, 1:2:end-1); dataOdor2 = dataPool(:, 2:2:end);

input.data{1, 1} = dataOdor1;  input.data{1, 2} = dataOdor2;

[trialNum, concNum] = size(dataOdor1);

resp = cell2mat(input.data);
dataY = [resp(:, 1:concNum); resp(:, concNum+1:2*concNum)];

conc = cell2mat(input.concList);
conc = conc';
dataX = [repmat(conc(1,:), [trialNum,1]); repmat(conc(2,:), [trialNum,1])];
dataX = log10(dataX);

% setup fitting function and method
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));
ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [-1 -1 -12]; % setup the range and initial value of the variable
opts.Upper = [10 15 1];
opts.StartPoint = [4 5 -6.5];

%% plot the curve on top of each other
figure; 

for i = 1 : trialNum
    plot(dataX(i, :), dataY(i, :), '*'); hold on;
end