function  sc = v2sc( v, dt, tBin, vmax )
%membran potential time corse to spike count

[N, numN] = size(v);

nBin = tBin/dt;
binCount = floor(N/nBin);

sc = zeros(binCount, numN);
for i = 1 : binCount
    tSection = (i-1)*nBin+1 : i*nBin;
    for j =1:numN
        ss = v(tSection, j);
        sc(i, j) = numel(find(abs(ss-vmax)<1E-8));
    end
end

end