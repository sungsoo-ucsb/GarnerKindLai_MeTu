function [ROI_mean, ROI_mean_temp, output_all, image_sum, total_volumes, vol_per_file] = analysis_by_stack_v4 (file_list_to_load, file_size, flcount, parameter_name, noise_threshold)

load(parameter_name)  
sliceNum = sliceNum*channel;
flyback = flyback*channel;
image_sum = zeros(pixel, pixel, channel, sliceNum);

%variables for image registration
usfac = 1;  %upsampling factor. Images will be registered to within 1/usfac of a pixel.
total_frames = sum(file_size(:,1))/channel;
vol_per_file = file_size(1,1)/(sliceNum+flyback);
total_volumes = floor(total_frames/(sliceNum+flyback))/channel;
usable_frames = total_frames - total_volumes*flyback;


disp('total frames:')
disp(total_frames)
disp('usable frames:')
disp(usable_frames)
disp('volumes per file:')
disp(vol_per_file)
disp('total volumes:')
disp(total_volumes)

%variables for ROI_signal_count
ROI_mean_temp = zeros(usable_frames, total_ROI, channel);
output_all = zeros (usable_frames, 4, channel);

import ScanImageTiffReader.ScanImageTiffReader;

count = 1;
count2 = 1;
for q = 1:flcount
    if q == flcount
        Tif_processed = zeros(pixel, pixel, floor(file_size(flcount,1)/(sliceNum+flyback))*sliceNum, channel);
    else
        Tif_processed = zeros(pixel, pixel, vol_per_file*sliceNum, channel);
    end
    file_name = file_list_to_load{q, 1};

    disp(['processing file ' num2str(q)])

    %step 1: discard the flyback frame(s)
    ch = 1;
    p = 1;
    vol = 1;
    A = sliceNum+flyback;
    index = linspace(1,A,A);
    z_plane_index = index;
    for i = 1:vol_per_file-1
        z_plane_index = [z_plane_index index];
    end 

    
    reader = ScanImageTiffReader(file_name);
    volume = reader.data();

    for i = 1:file_size(q,1)         
        if z_plane_index(i) < sliceNum + 1 % slice number
            image = volume(:,:,i);
            temp = reshape(volume(:,:,i), [1, pixel*pixel]);
            noise = max(mink(temp, floor(pixel*pixel*0.01)));

            for m = 1:pixel
                for n = 1:pixel
                    if image(m, n) > noise + noise_threshold
                        Tif_processed(n, m, p, ch) = image(m, n) - (noise_threshold + noise); 
                    end
                end
            end

            if channel == 1
                ch = 1;
                p = p+1;
            else
                if channel == 2 && ch == 1
                    ch = 2;
                else
                    ch = 1;
                    p = p+1;       
                end
            end
        end
        
        if rem(i,2000) == 0
            disp(['processed ' num2str(vol*2000) ' frames'])
            vol = vol +1;
        end    
    end

    if max(Tif_processed(:,:,end,channel), [], 'all') > 0
        disp(['length of Tif_processed: ' num2str(length(Tif_processed(pixel,pixel,:,channel)))])
    else
        disp('error')
    end

    %step 2
    Registered_tif = zeros (pixel, pixel, length(Tif_processed(pixel,pixel,:,channel)), channel);
    index = linspace(1,sliceNum,sliceNum);
    z_plane_index = index;
    for i = 1:vol_per_file-1
        z_plane_index = [z_plane_index index];
    end
    vol = 1;

    Tif_processed_gaussfilt = imgaussfilt3(Tif_processed(:,:,:,1), [0.1 0.1 sliceNum]); 
    

    for i = 1:length(Registered_tif(pixel,pixel,:,channel))
       
        for a = 1:channel
            if a == 1
                X = reference_image_name{z_plane_index(i),1};
                Fixed = imread(X);  %reference image
                if rgb(z_plane_index(i), 1) == 1
                    Fixed = rgb2gray(Fixed);
                end
            end
            if a == 2
                Fixed = imread('image_sum_a1_ch2.tif');  %reference image
                Fixed = rgb2gray(Fixed);
            end
            f = Fixed;

            Moving = Tif_processed_gaussfilt(:,:,i,1);

            g = Moving;
           
            if max(max(Moving)) > 0
                [output_temp Greg] = dftregistration(fft2(f),fft2(g),usfac);
                output_all(count, :, a) = output_temp; 
                row_shift = output_all(count,3,1);
                col_shift = output_all(count,4,1);
                Registered_tif(:,:,i,a) = circshift(Tif_processed(:,:,i,1), [row_shift col_shift]);
                count = count +1;
            else
                Registered_tif(:,:,i,a) = Tif_processed(:,:,i,1);
                count = count +1;
            end

        end

        if rem(i,2000) == 0
            disp(['processed ' num2str(vol*2000) ' frames'])
            vol = vol +1;
        end
    end

    if max(Registered_tif(:,:,end,channel), [], 'all') > 0
        disp(['length of Registered_tif: ' num2str(length(Registered_tif(pixel,pixel,:,channel)))])
    else
        disp('error in step 2')
    end
    
    clear Tif_processed
    clear Tif_processed_gaussfilt

    %step 3
    image_sum_tmp = zeros(pixel, pixel, channel, sliceNum);
    for a = 1:channel
        for i = 1:length(Registered_tif(pixel,pixel,:,channel))
            image(:,:) = Registered_tif(:,:,i,a);
            for m = 1:pixel
                for n = 1:pixel
                    image_sum_tmp(m,n,a,z_plane_index(i)) = image_sum_tmp(m,n,a,z_plane_index(i)) + image(m,n)/10;
                end
            end
        end
    end

    for a = 1:channel
        for i = 1:sliceNum
            image_sum(:,:,a,i) = image_sum(:,:,a,i) + image_sum_tmp(:,:,a,i);
        end
    end

    %step 4
    I = imread(X);
    for i = 1:length(Registered_tif(pixel,pixel,:,channel))
        for a = 1:channel
            image = Registered_tif(:,:,i,a);

            if i > 1
                image0 = Registered_tif(:,:,i-1,a);
            end
            if i <length(Registered_tif(pixel,pixel,:,channel))
                image1 = Registered_tif(:,:,i+1,a);
            end

            for j = 1:total_ROI
                mask = roipoly(I, x_cor(j,:), y_cor(j,:));
                temp = image(mask==1);
                if i >1
                    temp0 = image0(mask==1);
                end
                if i < length(Registered_tif(pixel,pixel,:,channel))
                    temp1 = image1(mask==1);
                end
                if ROI_inc_plane(j, z_plane_index(i)) == 1 && i>1 && i < length(Registered_tif(pixel,pixel,:,channel)) && z_plane_index(i) < 15  && z_plane_index(i) > 3 
                    ROI_mean_temp(count2,j,a) = mean((mean(temp) + mean(temp0) + mean(temp1)));  % use the average of 3 frames
                elseif ROI_inc_plane(j, z_plane_index(i)) == 1 && i>1 && i < length(Registered_tif(pixel,pixel,:,channel)) && z_plane_index(i) == 15
                    ROI_mean_temp(count2,j,a) = mean((mean(temp) + mean(temp0)));  % use the average of 2 frames
                elseif ROI_inc_plane(j, z_plane_index(i)) == 1 && i>1 && i < length(Registered_tif(pixel,pixel,:,channel)) && z_plane_index(i) == 2 
                    ROI_mean_temp(count2,j,a) = mean((mean(temp) + mean(temp1)));  % use the average of 2 frames
                elseif ROI_inc_plane(j, z_plane_index(i)) == 0 && count2 > 2
                    ROI_mean_temp(count2,j,a) = ROI_mean_temp(count2-1,j,a) ;
                elseif ROI_inc_plane(j, z_plane_index(i)) == 0 && count2 <= 2
                    ROI_mean_temp(count2,j,a) = 0;
                else
                    ROI_mean_temp(count2,j,a) = 0;
                end
            end
        end
        count2 = count2+1;
    end
    clear Registered_tif
end

%
ROI_mean = zeros(length(ROI_mean_temp), total_ROI, channel);
for i = 1:length(ROI_mean_temp)
    for a = 1:channel
        for j = 1:total_ROI
            ROI_mean(i, j, a) = ROI_mean_temp(i, j, a)-ROI_mean_temp(i, total_ROI, a); % minus the signal of the last ROI (background noise)
        end
    end
end


for i = 1:total_ROI-1
    ROI_mean(:,i,1) = movmean(ROI_mean(:,i,1), sliceNum, 1);
end


%================= find_peaks
I = imread(image_file_name);
ROI_mean_temp = zeros(total_ROI-1, sliceNum, channel);

ch=1;
for i = 1:sliceNum
    image = image_sum(:,:,ch,i);
    for j = 1:total_ROI
        mask = roipoly(I, x_cor(j,:), y_cor(j,:));
        temp = image(mask==1);
        ROI_mean_temp(j, i,ch) = mean(temp);
    end
end
