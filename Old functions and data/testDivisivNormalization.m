%% Backward ways:
%% PN uniform distribution
x = 0.01:0.01:1;
sigma = 0.05;

dHillEq = @(x) (2*x.^2 + sigma*(sigma*x).^0.5) ./ (x.^1.5 + sigma.^1.5).^2;

figure; 
plot(x, dHillEq(x)); 
xlabel('x'); ylabel('fraction'); title('Uniform')

%% PN exponential distribution
alpha = 1.25;
hillEq = @(x) (x.^1.5)./(x.^1.5 + sigma^1.5);

y = alpha * exp(-alpha * hillEq(x)) .* dHillEq(x);
figure; 
plot(x, y); 
xlabel('x'); ylabel('fraction'); title('Exp')

%% Forward way:
%% ORN exponential distribution
% generate 100 ORN activity data from exponential distribution


% define lateral inhibition
prop = 2.5/1000;
dNorm = @(x) (x.^1.5)./(x.^1.5 + sigma^1.5 + (prop*sum(x))^1.5);

%%
i = 0;
figure; set(gcf, 'Position', [680 10 800 760]);
for mu = [0.5 1 2 3]*sigma
    
    orns = random('Exponential', mu, [1000, 1]);
    pns = hillEq(orns);
    
    subplot(4, 4, i*4+1)
    histogram(orns, edges, 'Normalization', 'probability'); 
    xlabel('ORN'); ylabel('fraction'); title(['\lambda / \sigma= ', num2str(mu/sigma)]);

    subplot(4, 4, i*4+2)
    histogram(pns, edges, 'Normalization', 'probability'); 
    xlabel('PN-NoNorm'); ylabel('fraction'); 
    
    prop = 0.5/1000;
    dNorm = @(x) (x.^1.5)./(x.^1.5 + sigma^1.5 + (prop*sum(x))^1.5);
    pnnorms = dNorm(orns);
    subplot(4, 4, i*4+3)
    histogram(pnnorms, edges, 'Normalization', 'probability'); 
    xlabel('PN-DivNorm'); ylabel('fraction');  title(['m= ', num2str(prop)]);

    prop = 2.5/1000;
    dNorm = @(x) (x.^1.5)./(x.^1.5 + sigma^1.5 + (prop*sum(x))^1.5);
    pnnorms = dNorm(orns);
    subplot(4, 4, i*4+4)
    histogram(pnnorms, edges, 'Normalization', 'probability'); 
    xlabel('PN-DivNorm'); ylabel('fraction');  title(['m= ', num2str(prop)]);

    pause(0.1);i = i+1;
end
%% ORN distribution and PN distrinution using our results
lambda = 0.42;
sOffset = 9;
n = 1.4;

sens = random('Exponential', lambda, [1000, 1]);
sens = sens + sOffset;

edges = 0:0.02:1;

i = 0;
figure; set(gcf, 'Position', [680 10 800 760]);
for lc = -12 : 0.5 : -10.5
    ORNs = 1./(1+exp(-n*(lc + sens)));
    PNs = hillEq(ORNs);

    subplot(4, 4, i*4+1)
    histogram(ORNs, edges, 'Normalization','probability'); 
    xlabel('ORN'); ylabel('fraction'); title(['ln(c)=', num2str(lc)]); 

    subplot(4, 4, i*4+2)
    histogram(PNs, edges, 'Normalization','probability'); 
    xlabel('PN-NoNorm'); ylabel('fraction'); 

    prop = 0.5/1000;
    dNorm = @(x) (x.^1.5)./(x.^1.5 + sigma^1.5 + (prop*sum(x))^1.5);
    PNNorms = dNorm(ORNs);
    subplot(4, 4, i*4+3)
    histogram(PNNorms, edges, 'Normalization','probability'); 
    xlabel('PN-DivNorm'); ylabel('fraction'); title(['m= ', num2str(prop)]);
    
    prop = 2.5/1000;
    dNorm = @(x) (x.^1.5)./(x.^1.5 + sigma^1.5 + (prop*sum(x))^1.5);
    PNNorms = dNorm(ORNs);
    subplot(4, 4, i*4+4)
    histogram(PNNorms, edges, 'Normalization', 'probability'); 
    xlabel('PN-DivNorm'); ylabel('fraction');  title(['m= ', num2str(prop)]);
    
    pause(0.1); i = i+1;
end
