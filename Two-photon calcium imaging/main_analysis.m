clear

addpath(genpath('Z:\jennifer\Ring neurons\single_dot_blue\for github\Two-photon calcium imaging\Matlab function'))

temp_dir = dir('parameters_*.mat');
load(temp_dir.name)

parameter_name = temp_dir.name;


if isfile(['file_list_to_load_' dataset_name '.mat']) == 1
    load(['file_list_to_load_' dataset_name '.mat'])
else
    [file_list_to_load, file_size, flcount, file_offset] = import_tif_v4 (dataset_name, file_path, pixel, sliceNum, flyback);
end

tic

[ROI_mean, ROI_mean_temp, output_all, image_sum, total_volumes, vol_per_file] = analysis_by_stack_v4 (file_list_to_load, file_size, flcount, parameter_name, noise_threshold);

[F0, F_F0] = calculate_F_F0_v2 (ROI_mean, total_ROI);

[data_compressed, TP_count_compressed, TP_count_original, photodiode_count_compressed, stimuli_start_time, stim_time_method] = compress_data_file_v4(dataset_name, photodiode_threshold, total_trials, Dnum);

[data_compressed] = append_data_compressed_v9 (data_compressed, total_ROI, Dnum, F_F0, total_volumes, sliceNum, flyback);

[frame_weight, p_values] = statistical_test_by_stack_v10(Dnum, data_compressed, total_ROI, stimuli_start_time, dataset_name, image_sum, image_file_name, sliceNum, x_cor, y_cor);

save([file_path dataset_name '_analysis.mat'])

toc