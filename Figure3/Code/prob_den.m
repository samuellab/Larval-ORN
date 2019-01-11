function Ptrial = prob_den(rr,cc,AMean,Astd,x0Mean,x0Std,nMean,nStd,noiseStd)
    
M = size(rr,3);

%probabilities for generating latent variables x0, n given normal
%distributions for x0, n
Px0 = @(x0) 1./sqrt(2 * pi * x0Std.^2) * exp(-1/2*(x0 - x0Mean).^2/x0Std^2);
Pn = @(n) 1./sqrt(2 * pi * nStd.^2) * exp(-1/2*(n - nMean).^2/nStd^2);

%the probability for A is similar, but can be computed analytically...
%this mess below:
d = @(x0,n)1+exp(-n.*(cc - x0));
a = @(x0,n)1/(2*noiseStd^2) * sum(1./d(x0,n).^2,3) + 1/(2*Astd^2);
b = @(x0,n)1/noiseStd^2 * sum(rr./d(x0,n),3) + AMean/(Astd^2);
c = @(x0,n) -1/(2*noiseStd^2) * sum(rr.^2,3) - AMean^2/(2*Astd^2);
AIntegral = @(x0,n) 1/sqrt((2*pi*noiseStd^2)^M * 2*pi*Astd^2) *sqrt(pi./a(x0,n)).* ...
exp(b(x0,n).^2./(4.*a(x0,n))+c(x0,n));

%Probability density function for A, x0, n, with A integral
%precomputed analytically
PtrialFun = @(x0,n) AIntegral(x0,n) .* Px0(x0).* Pn(n);

abstol = 10^-6;
reltol = 10^-4;

%these functions are sharply peaked, so we need to choose numerical
%bounds carefully. We do this by finding the peak of the function,
%then integrating until the result stops increasing (it will do
%this because as the bound increase, the function will eventually fail to
%sample points on the peak of the function.
fun = @(parem)-PtrialFun(parem(1),parem(2));
[phat,pPeak] = fminsearch(fun,[x0Mean,nMean]);
widthx0 = x0Std;
pnew = 1;pold = 0;

%sample points around peak


%%

fac = 1;
stop = false;
while stop == false
pval = PtrialFun(phat(1)+fac*widthx0,phat(2));
if pval/(-pPeak) > 0.0001%increase the integral bounds
    fac = fac+1;
    stop = false;
    %pold = pnew;
else
    stop = true;
end
end
x0Lim = fac*widthx0;
%%


%%

%         fac = 1;
%         stop = false;
%         while stop == false
%             pnew = integral2(PtrialFun,phat(1)-fac*widthx0,phat(1)+fac*widthx0,...
%                 0,10,'AbsTol',abstol, 'RelTol',reltol);
%             if pnew > pold %increase the integral bounds
%                 fac = fac+1;
%                 stop = false;
%                 pold = pnew;
%             else
%                 stop = true;
%             end
%         end


%%
Ptrial = integral2(PtrialFun,phat(1)-x0Lim,phat(1)+x0Lim,...
    0,10,'AbsTol',abstol, 'RelTol',reltol);