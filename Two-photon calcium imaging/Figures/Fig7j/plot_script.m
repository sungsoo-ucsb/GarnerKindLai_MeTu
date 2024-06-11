clear
load('position_info.mat')
figure(1)
set(gcf,'color','w');
set(gcf, 'Position', [230,430,200,100])

for line = 1:2
    if line == 1 %R4d
        file_path_prefix = 'X:\Ring neurons\single_dot_blue\ER4d\';
        load('ER4d_list.mat')
        new_list = ER4d_list;
        load ellipse_data_ER4d.mat
    elseif line == 2 %R2
        file_path_prefix = 'X:\Ring neurons\single_dot_blue\ER2\';
        load('ER2_list.mat')
        new_list = ER2_list;
        load ellipse_data_ER2.mat
    end

subplot(1,2,line)
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
        outline1 = plot(x2, y2, 'color', line_color, 'LineWidth', 0.5, 'LineStyle','-');
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
        WXc = Wcentroid_x(1,k);
        WYc = Wcentroid_y(1,k);
        phi = deg2rad(orientation(1,k));

        t = linspace(0,2*pi,50);
        x = deg2rad(Xc) + deg2rad(a)*cos(t)*cos(phi) - deg2rad(b)*sin(t)*sin(phi);
        y = deg2rad(Yc) + deg2rad(a)*cos(t)*sin(phi) + deg2rad(b)*sin(t)*cos(phi);
        if line == 1
            MarkColor = [0 0.4 0.8];
        elseif line == 2
            MarkColor = [0.25 0.6 0];
        end

        ht = patch(rad2deg(x),rad2deg(y),MarkColor,'EdgeColor',MarkColor,'EdgeAlpha', 0.4,'FaceAlpha',0.1);
        %ht = patch(rad2deg(x),rad2deg(y),MarkColor,'EdgeColor','none','FaceAlpha',0.09);
        ht.Parent.LineWidth = 0.5;
        hold on

    end

    hold on
    axis equal
    ht.Parent.YLim = [-62 37];
    ht.Parent.XLim = [10 107];

    %
    x3 = [15 15]; y3 = [-60 40];
    y_axis_line = plot(x3, y3, '-k', 'LineWidth',0.5);

    for y_ticks = [  -60, -40, -20 , 0, 20, 40]
        x3 = [13 15]; y3 = [y_ticks y_ticks];
        y_tick_line = plot(x3, y3, '-k', 'LineWidth',0.5);
    end

    x3 = [20 100]; y3 = [-65 -65];
    x_axis_line = plot(x3, y3, '-k', 'LineWidth',0.5);

    for x_ticks = [20, 40, 60, 80, 100]
        x3 = [x_ticks x_ticks]; y3 = [-67 -65];
        x3_tick_line = plot(x3, y3, '-k', 'LineWidth',0.5);
    end
    axis off
    set(gca, 'box', 'off')

    ht.Parent.YLim = [-67 40];
    ht.Parent.XLim = [5 110];
    axis equal

    axis off
    set(gca, 'box', 'off')

end