function [A,AStd,n,nStd,noiseStd,functionVal,outputInfo] = ...
    FitMLE_global(r, conc, A0,AStd0,x0Vec0,x0StdVec0,n0,nStd0,noiseStd0)
% r (neuron response) should be a 3d array, rows should be concentration, columns should be
% trials, 3rd dimension is odor-ORN pairs. 
% AVec and other vectors are n-by-1 vector, n is the number of odor-neuron
% paris.

tic

%[~, ~, pairNum] = size(r);

fvalHistory = [];
% define the optimization function
pdfFun = @(parems) -likelihood_fun_global(parems(1), ...
    parems(2),x0Vec0, ...
    x0StdVec0, parems(3), ...
    parems(4), parems(5), r, conc);

%starting values for paremeters
paremStart = [A0; AStd0; n0; nStd0; noiseStd0];

options = optimset('TolFun',1e-2,'TolX',1e-2,'Display','none','MaxIter',100000,...
    'OutputFcn',@outfun,'MaxFunEvals',5000000);

[phat,functionVal,~,outputInfo] = fminsearch(pdfFun, paremStart, options);

A    = phat(1);
AStd = phat(2);
n = phat(3);
nStd = phat(4);
noiseStd = phat(5);

toc;

    function stop = outfun(x,optimValues,state)
%% display optimization progress
        stop = false;
        if strcmp(state,'iter') %&& mod(optimValues.iteration,10)==0
            
            fvalHistory = [fvalHistory; optimValues.fval];

            subplot(3,1,2)
            errorbar(optimValues.iteration,x(end-2),x(end-1),'bo')
            hold on;
            xlabel('iteration');ylabel('slope')

            subplot(3,1,3)
            errorbar(optimValues.iteration,x(1),x(2),'bo')
            hold on;
            xlabel('iter');ylabel('A');
            subplot(3,1,1)
            plot(optimValues.iteration,optimValues.fval,'bo')
            hold on;
            xlabel('iteration');ylabel('Function value')
            title(sprintf('Iteration %d',optimValues.iteration))
            drawnow;
%                 
%             
%             if fvalHistory(end) - fvalHistory(1) < -1e-4
%                 stop = true;
%                 
%             end
        
        end

    end
end