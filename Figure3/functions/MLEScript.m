% Warnning: run this code will take a very long time...
% final optimization result is saved as ./results/MLEFit.mat
clearvars;
[t123idx,sRow,sCol,concT123,rspT123] = prep_ds_data();

%% prepare the initial guess
load(which('initialParameters.mat'));
% 
pairNum = length(sRow);

for i = 1 : pairNum
    AVec0(i,1)  = aMatrix(sRow(i), sCol(i));
    x0Vec(i,1) = cMatrix(sRow(i), sCol(i));
end
% 
n = slop;

nStd = 1;
noiseStd = 1;

A = mean(AVec0);
AStd = std(AVec0);

x0StdVec = 1*ones(pairNum, 1);

%% fit
for iter = 1:5
    %first fit the global parameters
[A(iter+1),AStd(iter+1),n(iter+1),nStd(iter+1),noiseStd(iter+1),functionVal(iter+1),OutputInfo]...
    = FitMLE_global(rspT123, log10(concT123),...
    abs(normrnd(A(iter),AStd(iter)/2)),abs(normrnd(AStd(iter),AStd(iter)/2)),...
    x0Vec(:,iter),x0StdVec(:,iter),abs(normrnd(n(iter),nStd(iter))),abs(normrnd(nStd(iter),nStd(iter))),noiseStd(iter));
fprintf('Iter %d, A = %0.3f +- %0.3f, n = %0.3f +- %0.3f, function = %0.3f\n', iter, A(iter+1),...
    AStd(iter+1), n(iter+1),nStd(iter+1),functionVal(iter+1));
%then fit local parameters
[x0Vec(:,iter+1),x0StdVec(:,iter+1),funval(iter+1)] = ...
    FitMLE_local(rspT123,log10(concT123),A(iter+1),AStd(iter+1),x0Vec(:,iter),...
    x0StdVec(:,iter),n(iter+1),nStd(iter+1),noiseStd(iter+1));
fprintf('x0''s done, function = %0.3f\n',funval(iter+1))

end

%% Save fitted data 
% data is saved at ./results/MLEFit.mat.

%% generate 
