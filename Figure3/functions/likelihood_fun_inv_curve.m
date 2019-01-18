function [P] = likelihood_fun_inv_curve(Ai,x0i,ni,Amean,Astd,x0mean,x0std,nmean,nstd,noiseStd,r,conc)


er = @(A,x0,n) r - A./(1+exp(-n.*(conc - x0)));
Per = @(A,x0,n) sum(log(1./sqrt(2 * pi * noiseStd.^2))  - 1/2 * er(A,x0,n).^2./noiseStd.^2);
PA = @(Ai,Amean,Astd) log(1./sqrt(2 * pi * Astd.^2)) -1/2*(Ai - Amean).^2/Astd.^2;

P = Per(Ai,x0i,ni) + PA(Ai,Amean,Astd) + PA(x0i,x0mean,x0std) + PA(ni,nmean,nstd);
 



end