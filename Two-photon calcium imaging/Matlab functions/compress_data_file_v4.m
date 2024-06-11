function [data_compressed, TP_count_compressed, TP_count_original, photodiode_count_compressed, stimuli_start_time, stim_time_method] = compress_data_file_v4(dataset_name, photodiode_threshold, total_trials, Dnum)

load(['data_' dataset_name '.mat'])


timepoint = length(data(1,:));
data_compressed = zeros(Dnum, timepoint/10);
stim_time_method = 'photodiode';

TP_signal = zeros(1, length(data));
for i = 6:length(data)-10
    if  data(6,i) > 4 && data(6,i+1) > 4 && data(6,i+2) > 4
        TP_signal(1,i) = 5;
    else
        TP_signal(1,i) = 0;
    end
end


count = 1;
for i = 1:length(data_compressed(1,:))
    data_start_point = count;
    data_end_point = data_start_point+9;

    if i == 1
        data_compressed(1,i) = data(1, data_start_point);%time
        data_compressed(2,i) = data(2, data_start_point);
        data_compressed(3,i) = data(3, data_start_point);
        data_compressed(4,i) = data(4, data_start_point);
        data_compressed(5,i) = data(5, data_start_point);
        data_compressed(6,i) = 0;
        data_compressed(7,i) = 0;

    else
        data_compressed(1,i) = data(1, data_start_point);%time
        data_compressed(2,i) = mean(data(2, data_start_point:data_end_point));%LWB
        data_compressed(3,i) = mean(data(3, data_start_point:data_end_point));%L-R
        data_compressed(4,i) = mean(data(4, data_start_point:data_end_point));%Bhv_Cam_SYNC
        data_compressed(5,i) = mean(data(5, data_start_point:data_end_point));%new_trial
        if mean(TP_signal(1, data_start_point:data_end_point)) > 4 %TP_frame_SYNC
            data_compressed(6,i) = 5;
        else
            data_compressed(6,i) = 0;
        end
        data_compressed(7,i) = mean(data(7, data_start_point:data_end_point));%Photodiode
    end
    count = count+10;
end

photodiode_signal = zeros(1, length(data_compressed(1,:)));

if photodiode_threshold == -0.0606 || photodiode_threshold == -0.060|| photodiode_threshold == -0.046|| photodiode_threshold == -0.064
    photodiode_signal = movmean(data_compressed(7,:),4);
else
    if photodiode_threshold == -0.045 || photodiode_threshold == -0.05 || photodiode_threshold == -0.065 || photodiode_threshold == -0.07|| photodiode_threshold == -0.061
        photodiode_signal = movmean(data_compressed(7,:),6);
    else
        photodiode_signal = movmean(data_compressed(7,:),8);
    end
end

for i = 1:length(data_compressed(1,:))
    if photodiode_signal(1,i) < photodiode_threshold
        data_compressed(7,i) = 2;%Photodiode
    else
        if photodiode_signal(1,i) >= photodiode_threshold
            data_compressed(7,i) = 0;%Photodiode
        end
    end
end

TP_count_original = 0;
for i = 2:length(data(1,:))-1
    if TP_signal(1,i) > 4 && TP_signal(1,i-1) < 2
        TP_count_original = TP_count_original +1;
    end
end

TP_count_compressed = 0;
for i = 2:length(data_compressed(1,:))
    if data_compressed(6,i) > 4 && data_compressed(6,i-1) < 2
        TP_count_compressed = TP_count_compressed +1;
    end
end

for i = 1:length(data_compressed)
    if data_compressed(4,i) > 2.5
        data_compressed(4,i) = 3;
    else
         data_compressed(4,i) = 0;
    end
end


photodiode_count_compressed = 0;
for i = 2:length(data_compressed(1,:))
    if data_compressed(7,i) > 1 && data_compressed(7,i-1) < 1
        photodiode_count_compressed = photodiode_count_compressed +1;
    end
end

stimuli_start_time = zeros(total_trials, 1);

%find the starting time of the stimuli
j = 1;
for i = 5:length(data_compressed)-1
    if data_compressed(7, i) > 1  && data_compressed(7, i-1) < 1
        stimuli_start_time(j,1) = i;
        j = j+1;
        stim_time_method = 'photodiode';
    end
end
