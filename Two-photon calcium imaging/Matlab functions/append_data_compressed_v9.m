function [data_compressed] = append_data_compressed_v9 (data_compressed, total_ROI, Dnum, F_F0, total_volumes, sliceNum, flyback)

data_compressed(Dnum+1:Dnum+total_ROI, :) = zeros; %delta F_F0 of channel 1

A = sliceNum+flyback;
index = linspace(1,A,A);
z_plane_index = repmat(index, 1, total_volumes);


count = 1;
k = 1;
for i = 2:length(data_compressed)
    if count < length(z_plane_index)+1
        if data_compressed(6, i)<1 && data_compressed(6, i-1)>1
            if z_plane_index (count) < A-1 %TP_frame_SYNC
                for m = 1:total_ROI-1
                    data_compressed(Dnum+m, i) = F_F0(k, m);
                end
                count = count +1;
                k=k+1;
            else
                for m = 1:total_ROI-1
                    data_compressed(Dnum+m, i) = data_compressed(Dnum+m, i-1);
                end
                count = count +1;
            end
        else
            for m = 1:total_ROI-1
                data_compressed(Dnum+m, i) = data_compressed(Dnum+m, i-1);
            end
        end
    else
        break
    end

end

