function [x0fit,x0Std,funval] = FitMLE_local(r,concTotal,Astart,AstdStart,...
    x0start,x0stdStart,n,nstd,noiseStart)
% r should be a 3d array, rows should be concentration, columns should be odor/neuron pair, 
% 3rd dimension is trial. Assumes same number of trials for each odor
% neuron pair. 


r = permute(r,[1,3,2]);
concTotal = permute(concTotal,[1,3,2]);
h =waitbar(0,'Fitting individual curves...');
funval = 0;
for ii = 1:size(r,2)
%     tic%
    fvalHistory = [];
    
    rsmp = squeeze(r(:,ii,:));
    concentration = squeeze(concTotal(:,ii,:));
    tidx = isnan(rsmp(1,:));
    
    rsmp(:,tidx) = [];
    concentration(:,tidx) = [];
    pdfFun = @(parems) -likelihood_fun_local(parems(1),parems(2),Astart,AstdStart,...
        rsmp,n,nstd,concentration,noiseStart);%fit only x0 and noise
    
    %starting values for paremeters
    paremStart = [x0start(ii);x0stdStart(ii)];
    
    options = optimset('TolFun',1e-2,'Display','none','MaxIter',100000,'OutputFcn',@outfun,'MaxFunEvals',5000000);%,'PlotFcns',{@optimplotfval});   
    [phat,funcval,exitflag(ii)] = fminsearch(pdfFun,paremStart,options);
    x0fit(ii) = phat(1); x0Std(ii) = phat(2);

waitbar(ii/size(r,2),h);
    
funval= funval+funcval;
end
close(h);

%%
function stop = outfun(x,optimValues,state)

        stop = false;
        if strcmp(state,'iter') %&& mod(optimValues.iteration,10)==0
            
            fvalHistory = [fvalHistory; optimValues.fval];
            if fvalHistory(end) - fvalHistory(1) < -1e-4
                stop = true;
            end
        
        end

end

end