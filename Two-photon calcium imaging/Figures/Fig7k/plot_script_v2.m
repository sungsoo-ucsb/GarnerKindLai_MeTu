%%
clear
figure(1)
set(gcf,'color','w');
set(gcf, 'Position', [230,430,300,100])
ellipse_ratio = zeros(1,2);
ellipse_angle = zeros(1,2);
semimajor_axis = zeros(1,2);
semiminor_axis = zeros(1,2);
for line = 1:2
    subplot(1,2,line)
    if line == 1 %R4d
        load('ER4d_list.mat')
        load ellipse_data_ER4d.mat
        new_list = ER4d_list;
    elseif line == 2 %R2
        load('ER2_list.mat')
        load ellipse_data_ER2.mat
        new_list = ER2_list;
    end

    for i = 1:length(orientation)
        if orientation(1,i) < 0
            orientation(2,i)  = -(90-abs(orientation(1,i)));
        elseif orientation(1,i) > 0
            orientation(2,i)  = 90-(orientation(1,i));
        end
    end

    ax = gca;
    set(gca, 'box', 'off')
    for i = 1:length(orientation)
        x = orientation(2,i); y = major_axis_length(1,i)/minor_axis_length(1,i);
        ellipse_angle(i,line) = x; ellipse_ratio(i,line) = y;
        semimajor_axis(i,line) = major_axis_length(1,i); 
        semiminor_axis(i,line) = minor_axis_length(1,i); 
        dot1 = scatter(x,y,5,"filled");
        hold on
        if line == 1
            dot1.MarkerFaceColor = [0 0.4 0.8];
        elseif line == 2
            dot1.MarkerFaceColor = [0.25 0.6 0];
        end
    end
    dot1.Parent.XLim = [-95 95];
    dot1.Parent.YLim = [0.5 6.5];
    dot1.Parent.TickDir = 'out';
    dot1.Parent.XTick = [-90 -45 0 45 90];
    dot1.Parent.XTickLabel = [];
    dot1.Parent.YTick = [];
    dot1.Parent.YTickLabel = [];
    dot1.Parent.TickDir = 'out';
    dot1.Parent.LineWidth = 0.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.TickValues = [1 2 3 4 5 6];
end
