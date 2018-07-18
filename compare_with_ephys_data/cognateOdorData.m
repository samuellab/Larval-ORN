load('ORNSpikeData_Mathew13.mat')

t = ORNSpikeData;

[ornList, ~, ~] = unique(t(:, 1), 'stable');	%ORN name
conc = reshape(cell2mat(t(:, 3)), [5, length(ornList)]);
dff =  reshape(cell2mat(t(:, 4)), [5, length(ornList)]);

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
fprintf('%-5s\t%-5s\t%-5s\t%-5s\t\n', 'ORN', 'Amp', 'Slop', 'EC50');
for i = 1:19
    fprintf('%-5s\t%.2f\t%.2f\t%.2f\n', ...
        ornList{i}, ampVec, slop, kdVec(i));
end

dataXEn = conc' -  repmat(kdVec,  1, 5);
dataYEn = dff' ./ repmat(ampVec, 1, 5);

figure; 
plot(dataXEn', dataYEn', 'o'); hold on;
xPlot = linspace(min(dataXEn(:)), max(dataXEn(:)), 100);
yPlot = hillEq(1, slop, 0, xPlot);
plot(xPlot, yPlot, 'r'); xlabel('Relative Dose (log dilution)'); ylabel('Norm. firing rate');
hold off;

