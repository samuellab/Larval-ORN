function [Afit,AStd,x0fit,x0Std,noiseFit] = FitMLE(r,Astart,AstdStart,x0start,x0stdStart,noiseStart)
% r should be a 3d array, rows should be concentration, columns should be odor/neuron pair, 
% 3rd dimension is trial. Assumes same number of trials for each odor
% neuron pair. 

parfor ii = 1:size(r,2)
    tic
    rsmp = squeeze(r(:,ii,:));
    
    pdfFun = @(parems) -distrib_obj_fun_int(parems(1),parems(2),parems(3),parems(4),...
        rsmp,n,.2,CSample,parems(5),size(rsmp));
    
    %starting values for paremeters
    paremStart = [x0start(ii);x0stdStart;Astart;AstdStart;noiseStart];
    
    options = optimset('TolFun',1e-4,'Display','notify','MaxIter',10000,'MaxFunEvals',10000);%,'PlotFcns',@optimplotfval);
    [phat,~,exitflag(ii)] = fminsearch(pdfFun,paremStart,options);
    x0fit(ii) = phat(1); x0Var(ii) = phat(2)^2;
    Afit(ii) = phat(3); Avar(ii) = phat(4)^2;
    noiseFit(ii) = phat(5);
    toc;
    
end
x0Std = sqrt(x0Var);
AStd = sqrt(Avar);

%%
