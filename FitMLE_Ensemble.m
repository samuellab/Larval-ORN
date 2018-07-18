function [AVec,AStdVec,x0Vec,x0StdVec,n,nStd,noiseStd, fval] = FitMLE_Ensemble(r, conc, AVec0,AStdVec0,x0Vec0,x0StdVec0,n0,nStd0,noiseStd0)
% r should be a 3d array, rows should be concentration, columns should be
% trials, 3rd dimension is odor-ORN pairs. 
% AVec and other vectors are n-by-1 vector, n is the number of odor-neuron
% paris.

tic

[~, ~, pairNum] = size(r);

pdfFun = @(parems) -distrib_obj_fun_int_ensemble(parems(1 : pairNum), ...
    parems(pairNum+1 : pairNum*2), parems(2*pairNum+1 : pairNum*3), ...
    parems(3*pairNum+1 : pairNum*4), parems(pairNum*4 + 1), ...
    parems(pairNum*4 + 2), parems(pairNum*4 + 3), r, conc);

%starting values for paremeters
paremStart = [AVec0; AStdVec0; x0Vec0; x0StdVec0; n0; nStd0; noiseStd0];

options = optimset('TolFun',1e-4,'Display','notify','MaxIter',10000,'MaxFunEvals',10000  ,'PlotFcns',@optimplotfval);

[phat,fval ] = fminsearch(pdfFun, paremStart, options);

AVec    = phat(1 : pairNum);
AStdVec = phat(pairNum+1 : pairNum*2);
x0Vec   = phat(2*pairNum+1 : pairNum*3);
x0StdVec= phat(3*pairNum+1 : pairNum*4);
n = phat(pairNum*4 + 1);
nStd = phat(pairNum*4 + 2);
noiseStd = phat(pairNum*4 + 3);

toc;

