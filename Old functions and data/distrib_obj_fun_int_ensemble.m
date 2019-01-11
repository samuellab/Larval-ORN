function [P] = distrib_obj_fun_int_ensemble(AVec, AStdVec, x0Vec, x0StdVec, nMean, nStd, noiseStd, r, conc)

% r: concentratio number - by - trial number - by - pair number


% Ptrial = ones(1,size(data,2));
% C = permute(C,[3,2,1]);
Ppair = ones(size(r, 3),1);


parfor pp = 1:size(r, 3)

    
    cTrial = conc(:,:,pp);
    rTrial = r(:,:,pp);
    
    trialIdx = ~isnan(rTrial(1, :));
    
    cTrial = cTrial(:,trialIdx);
    rTrial = rTrial(: ,trialIdx);
    
    AMean = AVec(pp);       x0Mean = x0Vec(pp);
    Astd = AStdVec(pp);     x0Std = x0StdVec(pp);
    Ptrial = ones(size(rTrial, 2),1);
    for ii = 1:size(rTrial, 2) %loop over trials
        rr = permute(rTrial(:, ii),[3,2,1]); %one dose response curve for one trial
        cc = permute(cTrial(:, ii),[3,2,1]);

%         er = @(A,x0,n) rr ./ ( A./(1+exp(-n.*(cc - x0)))); % error function - difference between real data and estimated data using latent variables
        er = @(A,x0,n) rr  - A./(1+exp(-n.*(cc - x0))); % error function - difference between real data and estimated data using latent variables, additive noise

        %probabiliy of each trial given error fucntion with latent variable A's
        %and x0's. noiseStd represents trial-to-trial variability in response
        %for a single animal. 
%         PsAx0 = @(A,x0,n) prod(1./sqrt(2 * pi * noiseStd^2) * exp( - (er(A,x0,n)-1).^2./noiseStd^2),3); % muliplicative noise
        
%          PsAx0 = @(A,x0,n) prod(1./sqrt(2 * pi * (rr .* noiseStd).^2) .* exp( - er(A,x0,n).^2./(rr.*noiseStd).^2),3);% additive noise, propotional to the amp
         
          PsAx0 = @(A,x0,n) prod(1./sqrt(2 * pi * noiseStd^2) * exp( - er(A,x0,n).^2./noiseStd^2),3);% additive noise

        %probabilities for generating latent variables A, x0, given gaussian
        %model for A and x0
        PA = @(A) 1./sqrt( pi * Astd.^2) * exp(-(A - AMean).^2/Astd^2);
        Px0 = @(x0) 1./sqrt( pi * x0Std.^2) * exp(-(x0 - x0Mean).^2/x0Std^2);
        Pn = @(n) 1./sqrt( pi * nStd.^2) * exp(-(n - nMean).^2/nStd^2);
        
        PtrialFun = @(A,x0,n) PsAx0(A,x0,n) .* PA(A) .* Px0(x0).* Pn(n);
        Ptrial(ii) = integral3(PtrialFun,0,10,-10.5,-5,1,6,'RelTol',1e-2,'AbsTol',1e-4);
        
    end
    %%
    Ppair(pp) = sum(log(Ptrial));
    

end
P = sum(Ppair );