load(which('invMLEFitResults.mat'))
load(which('MLEFit.mat'))
%%
x0resid = cellfun(@(x0i,x0bar)x0i-x0bar,x0inv,num2cell(cMatrixMLE),'UniformOutput',false);
%%
Adist = cell2mat(cellfun(@(x)x(:),Ainv(:)','UniformOutput',false)');
ndist = cell2mat(cellfun(@(x)x(:),ninv(:)','UniformOutput',false)');
x0dist = cell2mat(cellfun(@(x)x(:),x0resid(:)','UniformOutput',false)');

%%
figure;histogram(Adist);
xlabel('A')
ylabel('Counts')
figure;histogram(ndist)
xlabel('n')
ylabel('Counts')
figure;histogram(x0dist)
xlabel('$x_0 - \bar{x_0}$','Interpreter','latex')
ylabel('Counts')