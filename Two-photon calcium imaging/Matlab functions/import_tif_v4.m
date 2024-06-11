function [file_list_to_load, file_size, flcount] = import_tif_v4 (dataset_name, file_path, pixel, sliceNum, flyback)

import ScanImageTiffReader.ScanImageTiffReader

file_size = [];
flcount = 0;

for i = 1:31
    file_name = [dataset_name  '_' num2str(i,'%05d') '.tif'];
    fn = fullfile(file_path, file_name);
    if exist(fn,'file')
        disp(['    reading a tiff file header: ' fn]);

        flcount = flcount+1;
        file_list_to_load{flcount,1} = fn;
        
        reader=ScanImageTiffReader(file_name);
        vol = reader.data();
        
        file_size(flcount,1) = length(vol(pixel,pixel,:));
        temp = reshape(vol(:,:,1), [1, pixel*pixel]);
        
    else
        break;
    end
end

A = sliceNum+flyback;
file_size(flcount,1) = floor(file_size(flcount,1)/A)*A;

save(['file_list_to_load_' dataset_name '.mat'],'file_list_to_load','file_size', 'flcount')

end

