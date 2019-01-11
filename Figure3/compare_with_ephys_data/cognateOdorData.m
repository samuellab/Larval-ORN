load('ORNSpikeData_Mathew13.mat')

t = ORNSpikeData;

[ornList, ~, ~] = unique(t(:, 1), 'stable');	%ORN name
conc = reshape(cell2mat(t(:, 3)), [5, length(ornList)]);
dff =  reshape(cell2mat(t(:, 4)), [5, length(ornList)]);

%% our data
rows = [3,  8,  6,  7,  9, 1, 16, 18, 17, 14, 10, 15, 5,  13, 11, 4   12];   
cols = [13, 16, 10, 14, 1, 4, 5,  17, 2,  8,  1,  6,  11, 20, 15, 12, 21];

lut = [1: 14, NaN, 15: 17, NaN];
% lut = [1, NaN, 2:7, NaN 9: 14, NaN, 15: 17, NaN];

load(fullfile('..', 'AnalysisResults', 'fitResults.mat'));

for i  = 1:length(rows)
    kdG(i) = cMatrix(rows(i), cols(i));   
end

%% pre fit

hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));

ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [0 0 -8]; % setup the range and initial value of the variable
opts.Upper = [250 10 0];
opts.StartPoint = [200 4 -6];

A0Vec = zeros(1, 19);  h0Vec = zeros(1, 19);	kd0 = zeros(1, 19);
fprintf('%-5s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

for i = 1:19
    
	xx = conc(:, i); 
	yy = dff(:, i);

    [fitresult, gof] = fit(xx, yy, ft, opts);   %fit

    rSq = gof.rsquare;  coeff = coeffvalues(fitresult);
    A0Vec(i) = coeff(1);  h0Vec(i)= coeff(2); kd0(i) = coeff(3);
    
    fprintf('%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
        ornList{i}, coeff(1), coeff(2), coeff(3), rSq);

end

%% 
h0 = median(h0Vec);  
A0 = median(A0Vec);
[slop, ampVec, kdVec, rSquare] = EnsembleMiniSearch(conc', dff', hillEq, h0, A0', kd0');

%%
disp('----------FIT CURVE ENSEMBLE:----------');
fprintf('%-5s\t%-5s\t\n', 'Slop', 'R^2');
fprintf('%.2f\t%.2f\t\n', slop, rSquare);

fprintf('Each individual curve:\n');
fprintf('%-5s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'ORN', 'Amp', 'Slop', 'EC50', 'EC50_GCaMP');

kdList = NaN(19, 1);
for i = 1:19
    if ~isnan(lut(i))
        kdList(i)= kdG(lut(i));
        fprintf('%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
            ornList{i}, ampVec, slop, kdVec(i), kdG(lut(i)));
    else
        fprintf('%-5s\t%.2f\t%.2f\t%.2f\n', ...
            ornList{i}, ampVec, slop, kdVec(i));
    end
end

dataXEn = conc' -  repmat(kdVec,  1, 5);
dataYEn = dff' ./ repmat(ampVec, 1, 5);

figure; 
plot(dataXEn', dataYEn', 'o'); hold on;
xPlot = linspace(min(dataXEn(:)), max(dataXEn(:)), 100);
yPlot = hillEq(1, slop, 0, xPlot);
plot(xPlot, yPlot, 'r'); xlabel('Relative Dose (log dilution)'); ylabel('Norm. firing rate');
hold off;


%%
figure; scatter(kdList(~isnan(kdList)), kdVec(~isnan(kdList))) ;
title(['corr.=', num2str(corr(kdList(~isnan(kdList)), kdVec(~isnan(kdList))))]);
