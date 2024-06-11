function [RF_weight, p_values] = statistical_test_by_stack_v10(Dnum, data_compressed, total_ROI, stimuli_start_time, dataset)

load(['data_' dataset '.mat'])

replicates = 10;
time_point = 2*1000*2;
positions = 38;
color = 1;
trials = replicates*positions*color;


RF_weight = zeros(positions, total_ROI);
p_values = zeros(positions, total_ROI);
for ROI = 1:total_ROI-1
    before_sti_mean = zeros(positions, trials);
    during_sti_mean = zeros(positions, trials);
    ROI_num = ROI; % deltaF_F0
    extracted_activity = zeros (trials, time_point);

    for i = 1:trials
        extracted_activity(i,:) = data_compressed(Dnum+ROI_num, (stimuli_start_time(i)-2*1000 : stimuli_start_time(i)-2*1000+time_point-1));
        before_sti_mean(1,i) =  mean(extracted_activity(i, 1001:2*1000));
        during_sti_mean(1,i) =  mean(extracted_activity(i, 2001:end));
    end


    for pos = 1:positions
        activity = zeros(replicates, 1);
        group1 = zeros(1,replicates);
        group2 = zeros(1,replicates);
        count = 1;
        for i = 1:trials
            if stimulus_info(3,i) == pos && stimulus_info(12,i) == 2
                activity(count,1) = during_sti_mean(1,i) - before_sti_mean(1,i);
                group1 (1,count) =  before_sti_mean(1,i);
                group2 (1,count) =  during_sti_mean(1,i);
                count = count + 1;
            end
        end
        temp = mean(activity(1:replicates, 1));
        RF_weight(pos,ROI) = temp; %blue
        p_values(pos,ROI) = signrank(group1(1,:), group2(1,:));
    end
end
