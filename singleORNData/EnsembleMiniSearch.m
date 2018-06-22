function [slop, ampVec, kdVec, rSquare] = EnsembleMiniSearch(dataX, dataY, hillEq, slop0, ampVec0, kdVec0)
% dataX, size: trialNum-by-xPoint, same as dataY.

[trialNum, ~] = size(dataX);
% slop0 = 4.4;
% % ampVec0 = 4 * ones(trialNum, 1);
% ampVec0 = [2.46; 2.49; 4.56; 5.57; 7.01; 5.39; 3.22; 4.99; 7.33; 6.45; 4.93; 5.07; ...
%     2.50; 2.50; 4.29; 5.68; 7.32; 5.62; 2.98; 4.24; 7.37; 6.50; 4.97; 4.86];
% % kdVec0  = -6 * ones(trialNum, 1);
% kdVec0 = [-5.96; -5.83; -5.69; -5.89; -5.82; -5.88; -5.78; -5.58; -5.71; -5.70; -5.95; ...
%     -6.10; -6.83; -6.76; -6.62; -6.78; -6.63; -6.81; -6.81; -6.53; -6.53; -6.56; ...
%     -6.75;-7.0];

p0 = [slop0; ampVec0; kdVec0];

eSSR = @(p) EnsembleSSR(p, dataX, dataY, hillEq);

% options = optimset('PlotFcns',@optimplotfval);
% [param, ssr]= fminsearch(eSSR, p0, options);

[param, ssr]= fminsearch(eSSR, p0);

slop= param(1);
ampVec = param(2:trialNum+1);
kdVec  = param(trialNum+2 : 2*trialNum+1);

sst = sum((dataY(:) - mean(dataY(:))).^2);
rSquare = 1-ssr/sst;

end

function ssr = EnsembleSSR(param, x, y, hillEq)
% input, output : 24*8
% x: parameters, one flop, 24 amp & 24 kd

[trialNum, samples] = size(x);

slop= param(1);
amp = param(2:trialNum+1);
kd  = param(trialNum+2 : 2*trialNum+1);

yHat = zeros(trialNum, samples);

for i = 1:trialNum
    yHat(i, :) = hillEq(amp(i), slop, kd(i), x(i, :));
end

ssr = sum((yHat(:) - y(:)).^2);

end