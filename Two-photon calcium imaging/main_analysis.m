
% Download a sample 2p dataset from https://sungsoo-nas1.mcdb.ucsb.edu:5001/sharing/nU6uifCJC
%(and all code in the "Two-photon calcium imaging" github folder)

%% Part 1: analysing 2p data
clear

addpath(genpath([pwd '\Matlab functions']))
addpath(genpath([pwd '\sample dataset']))

temp_dir = dir('sample dataset\parameters_*.mat');
load(temp_dir.name)

parameter_name = temp_dir.name;  
cd([pwd '\sample dataset'])
file_path = pwd;

tic
%read the tif files
if isfile(['file_list_to_load_' dataset_name '.mat']) == 1
    load(['file_list_to_load_' dataset_name '.mat'])
else
    [file_list_to_load, file_size, flcount] = import_tif_v4 (dataset_name, file_path, pixel, sliceNum, flyback);
end

%image registration, calculate F values of each ROI
[ROI_mean, ROI_mean_temp, output_all, image_sum, total_volumes, vol_per_file] = analysis_by_stack_v4 (file_list_to_load, file_size, flcount, parameter_name, noise_threshold);

%calculate dF/F
[F0, F_F0] = calculate_F_F0_v2 (ROI_mean, total_ROI);

%find the starting time of each visual stimuli (and compress the signal
%series for quicker processing)
[data_compressed, TP_count_compressed, TP_count_original, photodiode_count_compressed, stimuli_start_time, stim_time_method] = compress_data_file_v4(dataset_name, photodiode_threshold, total_trials, Dnum);

%synchronize dF/F with visual stimuli 
[data_compressed] = append_data_compressed_v9 (data_compressed, total_ROI, Dnum, F_F0, total_volumes, sliceNum, flyback);

%compare dF/F before and during the stimulations and conduct a statistical test
[RF_weight, p_values] = statistical_test_by_stack_v10(Dnum, data_compressed, total_ROI, stimuli_start_time, dataset_name);

save([file_path '\' dataset_name '_analysis.mat'])

toc

