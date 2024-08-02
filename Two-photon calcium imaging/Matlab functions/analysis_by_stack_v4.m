function [ROI_mean, ROI_mean_temp, output_all, image_sum, total_volumes, vol_per_file] = analysis_by_stack_v4 (file_list_to_load, file_size, flcount, parameter_name, noise_threshold)

load(parameter_name)  

image_sum = zeros(pixel, pixel, channel, sliceNum);

%variables for image registration
usfac = 1; 
total_frames = sum(file_size(:,1));
vol_per_file = file_size(1,1)/(sliceNum+flyback);
total_volumes = floor(total_frames/(sliceNum+flyback));
usable_frames = total_frames - total_volumes*flyback;

%variables
ROI_mean_temp = zeros(usable_frames, total_ROI);
output_all = zeros (usable_frames, 4);  %store all image registration outputs

import ScanImageTiffReader.ScanImageTiffReader.* 

count = 1;
count2 = 1;
for q = 1:flcount
    if q == flcount
        Tif_processed = zeros(pixel, pixel, floor(file_size(flcount,1)/(sliceNum+flyback))*sliceNum);
    else
        Tif_processed = zeros(pixel, pixel, vol_per_file*sliceNum);
    end
    file_name = file_list_to_load{q, 1};

    disp(['processing file ' num2str(q)])

    %step 1: discard the flyback frame(s)
    p = 1;
    vol = 1;
    A = sliceNum+flyback;
    index = linspace(1,A,A);
    z_plane_index = repmat(index, 1, vol_per_file);

    reader = ScanImageTiffReader(file_name);
    volume = reader.data();

    for i = 1:file_size(q,1)         
        if z_plane_index(i) < sliceNum + 1
            image = volume(:,:,i);
            temp = reshape(volume(:,:,i), [1, pixel*pixel]);
            noise = max(mink(temp, floor(pixel*pixel*0.01)));

            for m = 1:pixel
                for n = 1:pixel
                    if image(m, n) > noise + noise_threshold
                        Tif_processed(n, m, p) = image(m, n) - (noise_threshold + noise); 
                    end
                end
            end
            p = p+1;           
        end
        
        if rem(i,2000) == 0
            disp(['processed ' num2str(vol*2000) ' frames'])
            vol = vol +1;
        end    
    end

    if max(Tif_processed(:,:,end), [], 'all') > 0
        disp(['length of Tif_processed: ' num2str(length(Tif_processed(pixel,pixel,:)))])
    else
        disp('error')
    end

    %step 2
    Registered_tif = zeros(size(Tif_processed));
    index = linspace(1,sliceNum,sliceNum);
    z_plane_index = repmat(index, 1, vol_per_file);
    
    Tif_processed_gaussfilt = imgaussfilt3(Tif_processed(:,:,:), [0.1 0.1 sliceNum]); 
    vol = 1;
    for i = 1:length(Registered_tif(pixel,pixel,:))
        X = reference_image_name{z_plane_index(i),1};
        Fixed = imread(X);  %reference image
        if rgb(z_plane_index(i), 1) == 1
            Fixed = rgb2gray(Fixed);
        end

        f = Fixed;
        g = Tif_processed_gaussfilt(:,:,i,1); %moving image

        if max(max(g)) > 0
            [output_all(count, :)  ~] = dftregistration(fft2(f),fft2(g),usfac);
            row_shift = output_all(count,3);
            col_shift = output_all(count,4);
            Registered_tif(:,:,i) = circshift(Tif_processed(:,:,i), [row_shift col_shift]);
            count = count +1;
        else
            Registered_tif(:,:,i) = Tif_processed(:,:,i);
            count = count +1;
        end      

        if rem(i,2000) == 0
            disp(['processed ' num2str(vol*2000) ' frames'])
            vol = vol+1;
        end
    end

    if max(Registered_tif(:,:,end), [], 'all') > 0
        disp(['length of Registered_tif: ' num2str(length(Registered_tif(pixel,pixel,:)))])
    else
        disp('error in step 2')
    end
    
    clear Tif_processed
    clear Tif_processed_gaussfilt

    %step 3
    image_sum_tmp = zeros(pixel, pixel, sliceNum);

    for i = 1:length(Registered_tif(pixel,pixel,:))
        image(:,:) = Registered_tif(:,:,i);
        for m = 1:pixel
            for n = 1:pixel
                image_sum_tmp(m,n,z_plane_index(i)) = image_sum_tmp(m,n,z_plane_index(i)) + image(m,n)/10;
            end
        end
    end

    for i = 1:sliceNum
        image_sum(:,:,i) = image_sum(:,:,i) + image_sum_tmp(:,:,i);
    end

    %step 4
    I = imread(image_file_name);
    for i = 1:length(Registered_tif(pixel,pixel,:))   

        image = Registered_tif(:,:,i);
        if i > 1
            image0 = Registered_tif(:,:,i-1);
        end
        if i <length(Registered_tif(pixel,pixel,:))
            image1 = Registered_tif(:,:,i+1);
        end

        for j = 1:total_ROI
            mask = roipoly(I, x_cor(j,:), y_cor(j,:));
            if ROI_inc_plane(j, z_plane_index(i)) == 1
                temp = image(mask==1);
                if i > 1
                    temp0 = image0(mask==1);
                end
                if i < length(Registered_tif(pixel,pixel,:))
                    temp1 = image1(mask==1);
                end
            end

            if ROI_inc_plane(j, z_plane_index(i)) == 1 && z_plane_index(i) < 15  && z_plane_index(i) > 3
                ROI_mean_temp(count2,j) = mean((mean(temp) + mean(temp0) + mean(temp1)));  % use the average of 3 frames
            elseif ROI_inc_plane(j, z_plane_index(i)) == 1 && i>1  && z_plane_index(i) == 15
                ROI_mean_temp(count2,j) = mean((mean(temp) + mean(temp0)));  % use the average of 2 frames
            elseif ROI_inc_plane(j, z_plane_index(i)) == 1 && z_plane_index(i) == 2
                ROI_mean_temp(count2,j) = mean((mean(temp) + mean(temp1)));  % use the average of 2 frames
            elseif ROI_inc_plane(j, z_plane_index(i)) == 0 && count2 > 2
                ROI_mean_temp(count2,j) = ROI_mean_temp(count2-1,j) ;
            elseif ROI_inc_plane(j, z_plane_index(i)) == 0 && count2 <= 2
                ROI_mean_temp(count2,j) = 0;
            else
                ROI_mean_temp(count2,j) = 0;
            end
        end
        count2 = count2+1;
    end
    clear Registered_tif
end

%
ROI_mean = zeros(length(ROI_mean_temp), total_ROI);
for i = 1:length(ROI_mean_temp)
    for j = 1:total_ROI
        ROI_mean(i, j) = ROI_mean_temp(i, j)-ROI_mean_temp(i, total_ROI); % minus the signal of the last ROI (background noise)
    end
end

for i = 1:total_ROI-1
    ROI_mean(:,i) = sgolayfilt(ROI_mean(:,i),3,sliceNum);
end

