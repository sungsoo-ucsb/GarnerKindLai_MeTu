function [F0, F_F0] = calculate_F_F0_v2 (ROI_mean, total_ROI)

sort_list = zeros(length(ROI_mean), total_ROI);
q = 0.1; % lowest 10 percents
temp = round((q)*length(ROI_mean));
F0 = zeros(1, total_ROI);
F_F0 = zeros(length(ROI_mean), total_ROI);

for i = 1:total_ROI
    sort_list(:,i) = sort(ROI_mean(:,i), 'ascend');
end

for i = 1:total_ROI
    F0(1,i) = mean(sort_list(1:temp,i));
end

for i = 1:total_ROI
    F_F0(:,i) = (ROI_mean(:,i) - F0(1,i)) / F0(1,i); % deltaF / F0
end

