function [slop, ampVec, kdVec, rSquare] = EnsembleMiniSearch(dataX, dataY, hillEq, slop0, ampVec0, kdVec0)
% dataX, size: trialNum-by-xPoint, same as dataY.

[trialNum, ~] = size(dataX);

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