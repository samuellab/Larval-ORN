%% load data, averaged acorss trials of the raw data
% load(fullfile('.', 'data', 'AveRawDataMatrix.mat'));
load(fullfile('.', 'data', 'AveRawDataMatrix2ndRound.mat'));

%% plot the average activity vs concentration figure
[r, c, z] = size(dataMean);
odorConc = [1E-8 1E-7 1E-6 1E-5 1E-4];

% meanodorConc
actAveEachOdor = mean(dataMean, 2);
actAveEachOdor = reshape(actAveEachOdor, [r, z]);

actAveEachORN = mean(dataMean, 1);
actAveEachORN = reshape(actAveEachORN, [c, z]);

temp = reshape(dataMean, [r*c, z]);
actAveAllOdor = mean(temp, 1);

%standard error of the mean
actSEMEachOdor = std(dataMean,0, 2)./sqrt(c);
actSEMEachOdor = reshape(actSEMEachOdor, [r, z]);

actSEMEachORN = std(dataMean,0, 1)./sqrt(r);
actSEMEachORN = reshape(actSEMEachORN, [c, z]);

actSEMAllOdor = std(temp,0, 1)./sqrt(r*c);

%fit
Y = log10(actAveAllOdor');
X = ones(length(odorConc), 2);
X(:,2) = log10(odorConc');
slope = X\Y;

% plot
figure; 
errorbar(odorConc, actAveAllOdor, actSEMAllOdor, 'ok');
hold on

yy = X*slope;
plot(odorConc, 10.^yy, '--k');

title('Averaged ORN & Odor'); 
axis([5*10^-9 2*10^-4 0.03 1]);
xlabel('Concentration'); ylabel('\DeltaF/F');
set(gca,'XScale','log','YScale','log');
hold off

figure;
for i = 1:r
%     errorbar(odorConc, actAveEachOdor(i, :), actSEMEachOdor(i, :)); hold on;
    plot(odorConc, actAveEachOdor(i, :), 'o-'); hold on;
end
xlabel('Concentration'); ylabel('\DeltaF/F');
title('Each odor');  hold off; 
set(gca,'XScale','log','YScale','log'); 

figure;
for i = 1:c
    plot(odorConc, actAveEachORN(i, :), 'o-'); hold on;
end
xlabel('Concentration'); ylabel('\DeltaF/F');
title('Each ORN'); hold off;
set(gca,'XScale','log','YScale','log');