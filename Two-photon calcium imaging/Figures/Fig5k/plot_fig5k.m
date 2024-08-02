clear
load('position_info.mat')
figure(1)
set(gcf,'color','w');
set(gcf, 'Position', [230,430,200,400])

for line = 1:2
    if line == 1 %R4d
        file_path_prefix = 'X:\Ring neurons\single_dot_blue\ER4d\';
        load('ER4d_list.mat')
        new_list = ER4d_list; fly = 8;
        load ellipse_data_ER4d.mat
    elseif line == 2 %R2
        file_path_prefix = 'X:\Ring neurons\single_dot_blue\ER2\';
        load('ER2_list.mat')
        new_list = ER2_list; fly = 2;
        load ellipse_data_ER2.mat
    end

    subplot(2,1,line)

    vertex = [27	4.5
        27	-31.5
        45	-49.5
        81	-49.5
        99	-31.5
        99	4.5
        72	31.5
        45	22.5
        27	4.5
        ];
    for i = 1:8
        x2 = [vertex(i,1) vertex(i+1,1)]; y2 = [vertex(i,2) vertex(i+1,2)];
        line_color = [0 0 0 1];
        outline1 = plot(x2, y2, 'color', line_color, 'LineWidth', 1, 'LineStyle','-');
        %ht = patch(vertex(:,1),vertex(:,2), 'black','FaceAlpha',1); hold on
        hold on
    end

    for i = 1:length(orientation)
        if orientation(1,i) < 0
            orientation(2,i)  = -(90-abs(orientation(1,i)));
        elseif orientation(1,i) > 0
            orientation(2,i)  = 90-(orientation(1,i));
        end
    end

    for k = 1:length(major_axis_length)

        a = major_axis_length(1,k);
        b =  minor_axis_length(1,k);
        Xc = centroid_x(1,k);
        Yc = centroid_y(1,k);
        phi = deg2rad(orientation(1,k));

        t = linspace(0,2*pi,50);
        x = deg2rad(Xc) + deg2rad(a)*cos(t)*cos(phi) - deg2rad(b)*sin(t)*sin(phi);
        y = deg2rad(Yc) + deg2rad(a)*cos(t)*sin(phi) + deg2rad(b)*sin(t)*cos(phi);
        if line == 1
            if k == 19
                EdgeColor = 'red';
            elseif k == 20
                EdgeColor = 'green';
            elseif k == 21
                EdgeColor = 'blue';
            else
                EdgeColor = [0.4 0.4 0.4];
            end
             MarkColor = [0.6 0.6 0.6];
        elseif line == 2
            if k == 4
                EdgeColor = 'red';
            elseif k == 5
                EdgeColor = 'blue';
            else
                EdgeColor = [0.4 0.4 0.4];
            end
            MarkColor = [0.25 0.6 0]; 
        end

        ht = patch(rad2deg(x),rad2deg(y),MarkColor,'EdgeColor',EdgeColor,'EdgeAlpha', 0.7,'FaceAlpha',0,'LineWidth',1);
        ht.Parent.LineWidth = 1;
        hold on

    end

    hold on
    axis equal

    %
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
    axis off
    set(gca, 'box', 'off')

    ht.Parent.YLim = [-67 40];
    ht.Parent.XLim = [5 110];
    axis equal

    axis off
    set(gca, 'box', 'off')
end
%
% filename = 'fig5k';
% papersize = [2, 4];
% hf = gcf;
% units = 'inches';
% set(hf, 'PaperUnits', units)
% set(hf, 'papersize', papersize);
% set(hf, 'PaperPosition', [0, 0, papersize(1), papersize(2)]);
% print(hf, '-dpdf', filename);