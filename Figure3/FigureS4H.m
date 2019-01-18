load(fullfile('.', 'data', 'dataFigureS4H.mat'));

%% setup the fitting function
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));

ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [0 0 -8]; % setup the range and initial value of the variable
opts.Upper = [300 10 -1];
opts.StartPoint = [200 4 -6];

%% fit and plot
odor = 'ethyl acetate';
orn = {'Or42a', 'Or42b'};

fprintf('%25s\t%-5s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');
figure;

for i = 1:2
    if i == 1
        xx = Or42a(:, 1); 
        yy = Or42a(:, 2);
    elseif i==2
        xx = Or42b(:, 1); 
        yy = Or42b(:, 2);
    end
    [fitresult, gof] = fit(xx, yy, ft, opts);   %fit

    rSq = gof.rsquare; coeff = coeffvalues(fitresult);

    fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\t\n', 'ethyl acetate', ...
        orn{i}, coeff(1), coeff(2), coeff(3), rSq);

    xP = linspace(min(xx), max(xx), 50);
    yP = hillEq(coeff(1), coeff(2), coeff(3), xP);

    if i == 1
        plot(xx, yy, 'ok'); hold on;
    elseif i ==2
        plot(xx, yy, '.k', 'MarkerSize', 24); hold on;
    end
    
    plot(xP, yP, 'k'); xlabel('Dose(log dilution)'); ylabel('ORN response (spikes/s)');
end

title('dose-response to ethyl acetate');
legend({orn{1}, '', orn{2}}, 'Location', 'northwest');

% save figure
saveas(gcf, fullfile('.', 'results', 'figures', 'FigureS4H.fig'));