function [P] = likelihood_fun_local(x0Mean,x0Std,AMean,Astd,data,nMean,nStd,Conc,noiseStd)

%% data -- rows are concentrations, columns are trials. There's one x0 and
% one A guess for every trial

Ptrial = ones(1,size(data,2));

for ii = 1:size(data,2) %loop over trials
    rr = permute(data(:,ii),[3,2,1]); %one dose response curve for one trial
    cc = permute(Conc(:,ii),[3,2,1]);
    Ptrial(ii) = prob_den(rr,cc,AMean,Astd,x0Mean,x0Std,nMean,nStd,noiseStd);

end
%%

P = sum(log(Ptrial));

