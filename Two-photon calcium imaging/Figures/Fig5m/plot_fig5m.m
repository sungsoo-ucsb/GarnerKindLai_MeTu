clear
%
figure
set(gcf,'color','w');
set(gcf, 'Position', [230,430,200,400])

for line = 1:2
    subplot(2,1,line)
    if line == 1 %R4d
        load('X:\Ring neurons\single_dot_blue\analysis\ER4d_list')
        new_list = ER4d_list;
        load('ellipse_data_ER4d.mat','orientation')
    elseif line == 2 %R2
        load('X:\Ring neurons\single_dot_blue\analysis\ER2_list')
        new_list = ER2_list;
        load('ellipse_data_ER2.mat', 'orientation')
    end

    for i = 1:length(orientation)
        if orientation(1,i) < 0
            orientation(2,i)  = -(90-abs(orientation(1,i)));
        elseif orientation(1,i) > 0
            orientation(2,i)  = 90-(orientation(1,i));
        end
    end
    ellipse_angle= [];
    for i = 1:length(orientation)
        if orientation(2,i) < 0
            ellipse_angle(i) = 90+abs(orientation(2,i));
        elseif orientation(2,i) > 0
            ellipse_angle(i) = 90-abs(orientation(2,i));
        end
    end

    plot_angles = [0 20 40 60 80 100 120 140 160 180];
    plot_count = [];
    for i = 1:9
        count = 0;
        for j = 1:length(ellipse_angle)
            if ellipse_angle(j)>=plot_angles(i) &&  ellipse_angle(j)<plot_angles(i+1)
                count = count+1;
            end
        end
        plot_count(i) = count;
    end

    if line == 1
        plot_color = [0.7 0.7 0.7];
    elseif line == 2
        plot_color = [0.7 0.7 0.7];
    end
    
    for i = [0 30 60 90 120 150 180]
        theta =[]; rho =[];
        theta = [deg2rad(i) deg2rad(i)];
        rho(1,1) = 0;  rho(1,2) = 15;
        p = polarplot(theta,rho); rlim([0,15]);
        p.Color = [0 0 0];
        p.LineWidth =1;
        hold on
    end
    axis off
    if line == 1
        for i = [5 10 15]
            theta =[]; rho =[];
            theta = [linspace(0,pi,180)];
            rho(1,1:length(theta)) = i; 
            p = polarplot(theta,rho); rlim([0,15]);
            p.Color = [0 0 0];
            p.LineWidth = 1;
            hold on
        end
    elseif line == 2
        for i = [7/3 14/3 7]
            theta =[]; rho =[];
            theta = [linspace(0,pi,180)];
            rho(1,1:length(theta)) = i;
            p = polarplot(theta,rho); rlim([0,7]);
            p.Color = [0 0 0];
            p.LineWidth = 1;
            hold on
        end
    end
    for i = 1:9
        theta =[]; rho =[];
        theta = [deg2rad(plot_angles(i)) linspace(deg2rad(plot_angles(i)),deg2rad(plot_angles(i+1)),20)  deg2rad(plot_angles(i+1))];      
        rho = zeros(size(theta)) ;
        rho(1,1) = 0; rho(1, end) = 0; rho(1,2:end-1) = plot_count(i);
        p = polarplot(theta,rho); 
        if line == 1
            rlim([0,15]);
        elseif line == 2
            rlim([0,7]);
        end
        p.Color = plot_color;
        hold on
    end
    p.Parent.RGrid = 'off'; p.Parent.ThetaGrid = 'off'; 
    p.Parent.RTickLabel = [];
    p.Parent.ThetaTick = [0    30    60    90   120   150   180 ];
    p.Parent.ThetaTickLabel = [];
    p.Parent.ThetaLim = [0 360];

    ax_polar = gca;
    ax_cart = axes();
    ax_cart.Position = ax_polar.Position;
    for i = 1:9
        theta =[]; rho =[];
        theta = [deg2rad(plot_angles(i)) linspace(deg2rad(plot_angles(i)),deg2rad(plot_angles(i+1)),20)  deg2rad(plot_angles(i+1))];
        rho = zeros(size(theta)) ;
        rho(1,1) = 0; rho(1, end) = 0; rho(1,2:end-1) = plot_count(i);
        [xl,yl] = pol2cart(theta,rho);
        fill(xl,yl, plot_color)
        hold on
    end
    
    if line == 1
        xlim(ax_cart,[-15,15]);
        ylim(ax_cart,[-15,15]);
    elseif line == 2
        xlim(ax_cart,[-7,7]);
        ylim(ax_cart,[-7,7]);
    end
    axis square; set(ax_cart,'visible','off');
end
%
% filename = 'fig5m';
% papersize = [2, 4];
% hf = gcf;
% units = 'inches';
% set(hf, 'PaperUnits', units)
% set(hf, 'papersize', papersize);
% set(hf, 'PaperPosition', [0, 0, papersize(1), papersize(2)]);
% print(hf, '-dpdf', filename);