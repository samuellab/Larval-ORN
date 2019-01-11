function [slop, ampVec, kdVec, rSquare] = EnsembleMiniSearch4Cell(dataX, dataY, hillEq, slop0, ampVec0, kdVec0)
% both dataX and dataY is a trail-by-1 cell
% each elements of the cell contrain a vector, with different length

trialNum = length(dataX);

p0 = [slop0; ampVec0; kdVec0];

eSSR = @(p) EnsembleSSR(p, dataX, dataY, hillEq);

% options = optimset('PlotFcns',@optimplotfval);
% [param, ssr]= fminsearch(eSSR, p0, options);

[param, ssr]= fminsearch(eSSR, p0);

slop= param(1);
ampVec = param(2:trialNum+1);
kdVec  = param(trialNum+2 : 2*trialNum+1);

sst = sum(cellfun(@(x) sum((x-mean(x)).^2), dataY));

rSquare = 1-ssr/sst;
end

function ssr = EnsembleSSR(param, x, y, hillEq)
% input, output : 24*8
% x: parameters, one flop, 24 amp & 24 kd

trialNum = length(x);

slop= param(1);
amp = param(2:trialNum+1);
kd  = param(trialNum+2 : 2*trialNum+1);

yHat = cell(trialNum, 1);

ssr = 0;
for i = 1:trialNum
    yHat{i} = hillEq(amp(i), slop, kd(i), x{i});
    ssr = ssr + sum((yHat{i} - y{i}).^2);
end

end