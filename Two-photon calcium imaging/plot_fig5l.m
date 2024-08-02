clear

figure
set(gcf,'color','w');
set(gcf, 'Position', [230,430,200,400])
ellipse_ratio = zeros(1,2);
ellipse_angle = zeros(1,2);
semimajor_axis = zeros(1,2);
semiminor_axis = zeros(1,2);
for line = 1:2
    subplot(2,1,line)
    if line == 1 %R4d
        load('ER4d_list.mat')
        load ellipse_data_ER4d.mat
        new_list = ER4d_list;
    elseif line == 2 %R2
        load('ER2_list.mat')
        load ellipse_data_ER2.mat
        new_list = ER2_list;
        load ellipse_data_ER2.mat
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
        dot1 = scatter(x,y,10,"filled");
        hold on
        if line == 1
            if i == 19
                dot1.MarkerFaceColor = 'red';
            elseif i == 20
                dot1.MarkerFaceColor = 'green';
            elseif i == 21
                dot1.MarkerFaceColor = 'blue';
            else
                dot1.MarkerFaceColor = [0.3 0.3 0.3];
            end
        elseif line == 2
            if i == 4
                dot1.MarkerFaceColor = 'red';
            elseif i == 5
                dot1.MarkerFaceColor = 'blue';
            else
                dot1.MarkerFaceColor = [0.3 0.3 0.3];
            end
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
    dot1.Parent.LineWidth = 1;
    ax.YAxis.Color = 'k';
    ax.YAxis.TickValues = [1 2 3 4 5 6];
end
%
% filename = 'fig5l';
% papersize = [2, 4 ];
% hf = gcf;
% units = 'inches';
% set(hf, 'PaperUnits', units)
% set(hf, 'papersize', papersize);
% set(hf, 'PaperPosition', [0, 0, papersize(1), papersize(2)]);
% print(hf, '-dpdf', filename);
