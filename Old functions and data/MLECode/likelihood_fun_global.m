function [P] = likelihood_fun_global(AMean, Astd, x0Vec, x0StdVec, nMean, nStd, noiseStd, r, conc)

% r: concentratio number - by - trial number - by - pair number

Ppair = ones(size(r, 3),1);
%%

parfor pp = 1:size(r, 3)
    
    
    cTrial = conc(:,:,pp);
    rTrial = r(:,:,pp);
    
    trialIdx = find(~isnan(rTrial(1, :)));
    
    cTrial = cTrial(:,trialIdx);
    rTrial = rTrial(: ,trialIdx);
    
    x0Mean = x0Vec(pp);
    x0Std = x0StdVec(pp);
    Ptrial = ones(size(rTrial, 2),1);
    for ii = 1:size(rTrial, 2) %loop over trials for given odor-neuron pair
        rr = permute(rTrial(:, ii),[3,2,1]); %one dose response curve for one trial
        cc = permute(cTrial(:, ii),[3,2,1]);
        Ptrial(ii) = prob_den(rr,cc,AMean,Astd,x0Mean,x0Std,nMean,nStd,noiseStd);
    end
    %%
    Ppair(pp) = sum(log(Ptrial)); %sum the log of the probability
    
    
end
P = sum(Ppair );%sum the log of the probability