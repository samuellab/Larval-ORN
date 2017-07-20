%% load data, averaged acorss trials of the raw data
load(fullfile('.', 'data', 'AveRawDataMatrix.mat'));

%% plot the average activity vs concentration figure
[r, c, z] = size(dataRawAve);
odorConc = [1E-8 1E-7 1E-6 1E-5 1E-4];

% meanodorConc
actAveEachOdor = mean(dataRawAve, 2);
actAveEachOdor = reshape(actAveEachOdor, [r, z]);

actAveEachORN = mean(dataRawAve, 1);
actAveEachORN = reshape(actAveEachORN, [c, z]);

temp = reshape(dataRawAve, [r*c, z]);
actAveAllOdor = mean(temp, 1);

%standard error of the mean
actSEMEachOdor = std(dataRawAve,0, 2)./sqrt(c);
actSEMEachOdor = reshape(actSEMEachOdor, [r, z]);

actSEMEachORN = std(dataRawAve,0, 1)./sqrt(r);
actSEMEachORN = reshape(actSEMEachORN, [c, z]);

actSEMAllOdor = std(temp,0, 1)./sqrt(r*c);

% plot
figure; 
errorbar(odorConc, actAveAllOdor, actSEMAllOdor);
title('Averaged ORN & Odor'); 
xlabel('Concentration'); ylabel('\DeltaF/F');
set(gca,'XScale','log','YScale','log');

% fit the slope
Y = log10(actAveAllOdor');
X = ones(length(odorConc), 2);
X(:,2) = log10(odorConc');
slope = X\Y;


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