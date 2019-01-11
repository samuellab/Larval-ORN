function [P] = distrib_obj_fun_int(x0Mean,x0std,Amean,Astd,data,nmean,nstd,C,noiseStd,dataSz)

%% data -- rows are concentrations, columns are trials. There's one x0 and
% one A guess for every trial
%data = reshape(data,dataSz);
% Cmat = repmat(C,1,dataSz(2));
% Amat = repmat(As,dataSz(1),1);
% x0Mat = repmat(x0s,dataSz(1),1);
Ptrial = ones(1,size(data,2));
C = permute(C,[3,2,1]);

for ii = 1:size(data,2) %loop over trials
    d = permute(data(:,ii),[3,2,1]); %one dose response curve for one trial
    
    r = @(A,x0,n) d - A./(1+exp(-n.*(C - x0))); % error function - difference between real data and estimated data using latent variables
    
    %probabiliy of each trial given error fucntion with latent variable A's
    %and x0's. noiseStd represents trial-to-trial variability in response
    %for a single animal. 
    PsAx0 = @(A,x0,n) prod(1./sqrt(2 * pi *noiseStd^2) * exp( - r(A,x0,n).^2./noiseStd^2),3);

    %probabilities for generating latent variables A, x0, given gaussian
    %model for A and x0
    PA = @(A) 1./sqrt( pi * Astd.^2) * exp(-(A - Amean).^2/Astd^2);
    Px0 = @(x0) 1./sqrt( pi * x0std.^2) * exp(-(x0 - x0Mean).^2/x0std^2);
    Pn = @(n) 1./sqrt( pi * nstd.^2) * exp(-(n - nmean).^2/nstd^2);
    PtrialFun = @(A,x0,n) PsAx0(A,x0,n) .* PA(A) .* Px0(x0).* Pn(n);
    Ptrial(ii) = integral3(PtrialFun,0,10,-10,-2,1,8,'RelTol',1e-2,'AbsTol',1e-4);
    
    
end
%%
P = log(prod(Ptrial));