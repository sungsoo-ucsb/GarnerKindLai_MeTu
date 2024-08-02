clear
load('20240207f_00005_analysis.mat', 'Dnum', 'data_compressed','stimuli_start_time','total_trials')
load('ROI_cor_20240207f_00005.mat')
load('data_20240207f_00005.mat', 'stimulus_info')
load('position_info.mat')
%
figure
set(gcf,'color','w');
set(gcf, 'Position', [230,430,350,100])
set(gca, 'box', 'off')

for ROI = 1:3

    replicates = 10;
    extracted_activity = zeros(replicates, 6000);
    time_point = 6000;
    
    count = 1;
    for i = 1:total_trials
        if ROI ==1
            pos = 5;
        elseif ROI == 2
            pos = 36;
        elseif ROI == 3
            pos = 20 ;
        end
        
        if stimulus_info(3,i) == pos
            extracted_activity(count,:) = data_compressed(Dnum+ROI, (stimuli_start_time(i)-2*1000 : stimuli_start_time(i)-2*1000+time_point-1));
            count = count +1;
        end    
    end

    %stat test
    before_stim = zeros(1,replicates);
    during_stim = zeros(1,replicates);
    diff = zeros(1,replicates);
    for i = 1:10
        before_stim (1,i) = mean(extracted_activity(i,1001:2000));
        during_stim (1,i) = mean(extracted_activity(i,2001:4000));
        diff(1,i) =  during_stim(1,i) - before_stim(1,i);
    end
    p_value1 = signrank(before_stim(1,:), during_stim(1,:));

    if ROI == 1
        subplot(1,3,1)
    elseif ROI == 2
        subplot(1,3,3)
    elseif ROI == 3
        subplot(1,3,2)
    end
 
    ax = gca; ax.LineWidth = 1;
    hold on
    x1 = []; x2 = [];
    x1(1:10) = 1; x2(1:10) = 2;
    y1 = []; y2 = [];
    for i = 1:10
        y1(i) = before_stim(1,i);
        y2(i) = during_stim(1,i);
    end
    sz = 14;

    %add a gray line connecting the dots
    for i = 1:length(y2)
        plot([x1(i) x2(i)], [y1(i) y2(i)], 'Color', [0.8 0.8 0.8])
    end
    
    if ROI == 1
        MarkColor = "red";
    elseif ROI == 2
        MarkColor = "green";
    elseif ROI == 3
        MarkColor = "blue";
    end
        t1 = scatter(x1,y1,sz,MarkColor,'filled','MarkerFaceAlpha', 0.5, 'jitter','on','jitterAmount',0.15, 'MarkerEdgeColor', MarkColor);
        t2 = scatter(x2,y2,sz,MarkColor,'filled','MarkerFaceAlpha', 0.5, 'jitter','on','jitterAmount',0.15, 'MarkerEdgeColor', MarkColor);
    
    
    if ROI == 1 || ROI == 3
        YLim_up = 2.7;
        t2.Parent.YLim = [-0.2 YLim_up]; t2.Parent.XLim = [0 2.5];
    elseif ROI == 2
        YLim_up = 1.5;
        t2.Parent.YLim = [-0.2 YLim_up]; t2.Parent.XLim = [0 2.5];
    end

    x = [1  2]; y = [YLim_up YLim_up];
    xline = plot(x, y, '-k', LineWidth=1);
    x = [0.5 0.5]; y = [-0.2 YLim_up];
    xline = plot(x, y, '-k', LineWidth=1);
    x = [0.4 0.5]; y = [0 0];
    yline = plot(x, y, '-k', LineWidth=1);
    x = [0.4 0.5]; y = [1 1];
    yline = plot(x, y, '-k', LineWidth=1);
    x = [0.4 0.5]; y = [2 2];
    yline = plot(x, y, '-k', LineWidth=1);

    hold on
    axis off
end

% filename = 'fig5g';
% papersize = [3.5, 1];
% hf = gcf;
% units = 'inches';
% set(hf, 'PaperUnits', units)
% set(hf, 'papersize', papersize);
% set(hf, 'PaperPosition', [0, 0, papersize(1), papersize(2)]);
% print(hf, '-dpdf', filename);
