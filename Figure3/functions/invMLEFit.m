
clearvars
load(which('MLEFit.mat'));
[t123idx,sRow,sCol,concT123,rspT123] = prep_ds_data();


x0Vec = cMatrixMLE(t123idx);

%% perform fitting on individual curves using previously generated distributions
h = waitbar(0,sprintf('Fitting odor-ORN pair %d/%d',0,size(rspT123,3)));
for ii = 1:size(rspT123,3)
    data = rspT123(:,:,ii);
    conc = log10(concT123(:,:,ii));
    data(:,isnan(data(1,:)))= [];
    conc(:,isnan(conc(1,:))) = [];
    for t = 1:size(data,2)
        optfun = @(parems) -likelihood_fun_inv_curve(parems(1),parems(2),parems(3),...
            A,AStd,x0Vec(ii),x0Std,n,nStd,noiseStd,...
            data(:,t),conc(:,t));
        
        paremStart = [A,x0Vec(ii),n];
        
        options = optimset('TolFun',1e-6,'Display','notify');
        [phat,funval,exitflag{ii}(t)] = fminsearch(optfun,paremStart,options);

        Ai{ii}(t) = phat(1);
        x0i{ii}(t) = phat(2);
        ni{ii}(t) = phat(3);
%         
    end
    waitbar(ii/size(rspT123,3),h,sprintf('Fitting odor-ORN pair %d/%d',ii,size(rspT123,3)));
end
close(h);


%% reformat results 

Ainv = cell(size(t123idx));
Ainv(t123idx) = Ai;

ninv = cell(size(t123idx));
ninv(t123idx) = ni;

x0inv = cell(size(t123idx));
x0inv(t123idx) = x0i;
saveFile = which('invMLEFitResults.mat');
save(saveFile,'Ainv','ninv','x0inv')