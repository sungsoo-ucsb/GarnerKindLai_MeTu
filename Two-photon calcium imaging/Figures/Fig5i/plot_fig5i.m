clear
load('position_info.mat')
%
figure
set(gcf,'color','w');
set(gcf, 'Position', [230,430,200,400])

for line = 1:2
    if line == 1
        load('20240207f_00005_analysis.mat', 'Dnum', 'data_compressed','stimuli_start_time','total_trials','RF_weight','p_values','total_ROI')
        load('data_20240207f_00005.mat', 'stimulus_info')
        load('ER4d_list'); new_list = ER4d_list;
        fly = 8;
    elseif line == 2
        load('20240204c_00002_analysis.mat', 'Dnum', 'data_compressed','stimuli_start_time','total_trials','RF_weight','p_values','total_ROI')
        load('data_20240204c_00002.mat', 'stimulus_info')
        load('ER2_list'); new_list = ER2_list;
        fly = 2;
    end
    subplot(2,1,line)

    x = position_info(:,1); y = position_info(:,2);
    v = zeros(size(RF_weight));
    for ROI = 1:total_ROI-1
        for pos = 1:38
            if p_values(pos,ROI) < 0.05 && RF_weight(pos,ROI) > max(RF_weight(:,total_ROI-1)) && RF_weight(pos,ROI) > 0.2*max(RF_weight(:,ROI))
                v(pos, ROI) = RF_weight (pos,ROI);
            end
        end
    end
    grid_x = repmat(linspace(27,99,55), 55, 1);
    grid_y = repmat(linspace(31.5, -49.5,55)', 1,55);
    grid_z = zeros(55,55);

    for ROI = 1:length(new_list{fly,2})
        grid_z = zeros(55,55);
        for k = 1:38
            for i = 1:length(grid_x)
                for j = 1:length(grid_y)
                    if grid_x(i,j) == x(k,1) && grid_y(i,j) == y(k,1)
                        grid_z(i, j) = v(k,ROI);
                    end
                end
            end
        end

        for i = 1:length(grid_x)
            for j = 1:length(grid_y)
                if grid_z(i,j) == 0
                    grid_z(i,j) = griddata(x,y,v(:,ROI), grid_x(i,j), grid_y(i,j));
                end
            end
        end

        [ct, ht] = contour(grid_x, grid_y, grid_z);
        ht.LevelList = [0.2*max(grid_z, [],"all")];
        ht.EdgeColor  = 'none';
        if line == 1
            if ROI == 1
                ht.FaceColor = 'red';
            elseif ROI == 2
                ht.FaceColor = 'green';
            elseif ROI == 3
                ht.FaceColor = 'blue';
            end
        elseif line == 2
            if ROI == 1
                ht.FaceColor = 'red';
            elseif ROI == 2
                ht.FaceColor = 'blue';
            elseif ROI == 3
                ht.FaceColor = 'blue';
            end
        end
        ht.FaceAlpha = 0.3;
        hold on
        axis equal
        [ct, ht] = contour(grid_x, grid_y, grid_z);
        ht.LevelList = [0.5*max(grid_z, [],"all")];
        ht.EdgeColor  = 'none';
        if line == 1
            if ROI == 1
                ht.FaceColor = 'red';
            elseif ROI == 2
                ht.FaceColor = 'green';
            elseif ROI == 3
                ht.FaceColor = 'blue';
            end
        elseif line == 2
            if ROI == 1
                ht.FaceColor = 'red';
            elseif ROI == 2
                ht.FaceColor = 'blue';
            elseif ROI == 3
                ht.FaceColor = 'blue';
            end
        end
        ht.FaceAlpha = 0.4;
        hold on

        if ROI == 1
            ht.Parent.YLim = [-67 40];
            ht.Parent.XLim = [5 110];
            x3 = [15 15]; y3 = [-60 40];
            y_axis_line = plot(x3, y3, '-k', 'LineWidth',1);

            for y_ticks = [  -60, -40, -20 , 0, 20, 40]
                x3 = [13 15]; y3 = [y_ticks y_ticks];
                y_tick_line = plot(x3, y3, '-k', 'LineWidth',1);
            end

            x3 = [20 100]; y3 = [-65 -65];
            x_axis_line = plot(x3, y3, '-k', 'LineWidth',1);

            for x_ticks = [20, 40, 60, 80, 100]
                x3 = [x_ticks x_ticks]; y3 = [-67 -65];
                x3_tick_line = plot(x3, y3, '-k', 'LineWidth',1);
            end
        end
        axis off
        ht.Parent.XTickLabel = '';
        ht.Parent.YTickLabel = '';
        ht.Parent.TickDir = 'out';
        ht.Parent.LineWidth = 1;
        set(gca, 'box', 'off')
    end
end
%
% filename = 'fig5i';
% papersize = [2, 4];
% hf = gcf;
% units = 'inches';
% set(hf, 'PaperUnits', units)
% set(hf, 'papersize', papersize);
% set(hf, 'PaperPosition', [0, 0, papersize(1), papersize(2)]);
% print(hf, '-dpdf', filename);