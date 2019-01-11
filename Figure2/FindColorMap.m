function [colorMapSort] = FindColorMap()

load(fullfile('..', 'Figure4','data', 'log10EC50.mat'));
% cMatrix(23, 12) = -9;   %2-heptanone, 85c
% cMatrix(24, 16) = -9.5;  %methyl salicylate, 22c 
log10EC50(isnan(log10EC50)) = 0;
[~, score, ~,~,~, ~] = pca(-log10EC50);
[~, I] = sort(score(:, 1));

myColorMap = jet(34);
colorMapSort = zeros(34, 3);
for i =1:34
    colorMapSort(I(i), :) = myColorMap(i, :);
end

end

%% 
% load(dataFileName );
% data = dataMean;
%  
% 
% %% PCA
% dataTall = permute(data, [1 3 2]); 
% dataTall = reshape(dataTall,[],size(data,2),1);
% dataTall(171, :) = 0;
% [~, score, ~, ~, ~, ~] = pca (dataTall);
% 
% %% 
% points = score(5*(1:34), 1:3);
% center = score(171, 1:3);
% points = points - repmat(center, [34, 1]);
%  
% r_theta_phi = zeros(34, 3);
% for i = 1:34
%     x = points(i,1);
%     y = points(i,2);
%     z = points(i,3);
%     r_theta_phi(i,1) = sqrt(x^2 + y^2 + z^2); %radius, r
%     r_theta_phi(i,2) = acos(z/sqrt(x^2 + y^2 + z^2)); %inclination, theta
%     r_theta_phi(i,3) = atan(y/x); %azimuth, phi
% end
% 
% %% calculate the angle
% directionCenter = [mean(r_theta_phi(:,2))  mean(r_theta_phi(:,3))];
% direction = [r_theta_phi(:,2)-directionCenter(1) r_theta_phi(:,3)-directionCenter(2)];
% 
% rot = zeros(1,34);
% for i =1:34
%     rot(i) = atan2(direction(i, 2), direction(i, 1));
% end
% 
% [~, I] = sort(rot);
% 
% %%
% myColorMap = hsv(34);
% colorMapSort = zeros(34, 3);
% for i =1:34
%     colorMapSort(I(i), :) = myColorMap(i, :);
% end
% 
% %%
% figure;
% for i = 1:34
%      plot(r_theta_phi(i,2), r_theta_phi(i,3), 'o', 'MarkerFaceColor', colorMapSort(i, :));
%      hold on
% end
% plot(mean(r_theta_phi(:,2)), mean(r_theta_phi(:,3)), 'r+');
% 
% %add odor name
% dx = 0.01; dy = 0.01;  % displacement so the text does not overlay the data points
% text(r_theta_phi(:,2) + dx, r_theta_phi(:,3) +dy, odorList);  
% hold off;
% end