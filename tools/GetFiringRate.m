function [ tranM, aveFR ] = SimplifyResults( sc, tBin, pDuration )
% reSimplify the results, detect if the response transient, and get the
% averaged firing rate. 

% tBin = 0.1
% pDuration = 0.5

[tSeg, numN] = size(sc);
tPoint = tSeg * tBin / pDuration;
range = pDuration/tBin;

tranM = zeros(tPoint, numN);
aveFR = zeros(tPoint, numN);

for i = 1 : tPoint
    for n = 1 : numN
        iL = (i-1)*range + 1;
        iR = i*range;
        sig = sc(iL:iR, n);
        if sig(1)>sum(sig(2:end))
            tranM(i, n) = 1;
            aveFR(i, n) = sig(1)/tBin;
        else
            tranM(i, n) = 0;
            aveFR(i, n) = sum(sig)/pDuration;
        end
    end
end

end