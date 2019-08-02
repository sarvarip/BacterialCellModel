close all
clear
load('FYP_04_06_par_ex_15')
idx = 1;
aac_array = [7500, 300, 300, 300];
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:10:4000
    ax1 = axes('Position',[0.3 0.15 0.7 0.7],'Box','off');
    ax2 = axes('Position',[0.1 0.7 0.25 0.2],'Box','off');
    ax3 = axes('Position',[0.1 0.1 0.25 0.2],'Box','off');
    ax = axes('Position',[0 0 1 1],'Box','off');
    if P_count_vec_array(i, 4) == 0
        axes(ax1)
        label = {'R-proteins: ', 'Enzyme: ', 'Housekeeping: '};
        p= pie(P_count_vec_array(i, :).*aac_array);
        t = p(2);
        t.FontSize = 15;
        t = p(4);
        t.FontSize = 15;
        t = p(6);
        t.FontSize = 15;
        colormap([1 0 0;      %// red
              0 1 0;      %// green
              0 0 1]);      %// blue
        hText = findobj(p,'Type','text'); % text object handles
        percentValues = get(hText,'String'); % percent values
        combinedtxt = strcat(label,percentValues'); % strings and percent values
        hText(1).String = combinedtxt(1);
        hText(2).String = combinedtxt(2);
        hText(3).String = combinedtxt(3);
    else
        axes(ax1)
        label = {'R-proteins: ', 'Enzyme: ', 'Housekeeping: ', 'Heterologous: '};
        p= pie(P_count_vec_array(i, :).*aac_array);
        t = p(2);
        t.FontSize = 15;
        t = p(4);
        t.FontSize = 15;
        t = p(6);
        t.FontSize = 15;
        hText = findobj(p,'Type','text'); % text object handles
        percentValues = get(hText,'String'); % percent values
        combinedtxt = strcat(label,percentValues'); % strings and percent values
        hText(1).String = combinedtxt(1);
        hText(2).String = combinedtxt(2);
        hText(3).String = combinedtxt(3);
        hText(4).String = combinedtxt(4);
        t = p(8);
        t.FontSize = 15;
        colormap([1 0 0;      %// red
          0 1 0;      %// green
          0 0 1;      %// blue
          .5 .5 .5])  %// grey
    end
    title(['Proteome partition, time = ', num2str(time_ss(i)), ' mins'], 'FontSize', 15);
    axes(ax2);
    scatter(time_ss(1:10:i), energy_array(1:10:i), 'c')
    ylim([0 8*10^5])
    title('Energy', 'FontSize', 15);
    ylabel('a.u.', 'FontSize', 15);
    xlabel('Time (mins)', 'FontSize', 15);
    axes(ax3);
    scatter(time_ss(1:10:i), total_inst_gr_array(1:10:i), 'm')
    title('Inst. Growth Rate', 'FontSize', 15);
    ylabel('Hr^-^1', 'FontSize', 15);
    xlabel('Time (mins)', 'FontSize', 15);
    ylim([0 3])
    fig = gcf;
    fig.CurrentAxes = ax;
    M(idx) = getframe;
    ax.Units = 'normalized';
    disp(i)
    idx = idx+1;
    %close
end
for i = 4100:100:10000
    ax1 = axes('Position',[0.3 0.15 0.7 0.7],'Box','off');
    ax2 = axes('Position',[0.1 0.7 0.25 0.2],'Box','off');
    ax3 = axes('Position',[0.1 0.1 0.25 0.2],'Box','off');
    ax = axes('Position',[0 0 1 1],'Box','off');
    axes(ax1)
    label = {'R-proteins: ', 'Enzyme: ', 'Housekeeping: ', 'Heterologous: '};
    p= pie(P_count_vec_array(i, :).*aac_array);
    t = p(2);
    t.FontSize = 15;
    t = p(4);
    t.FontSize = 15;
    t = p(6);
    t.FontSize = 15;
    t = p(8);
    t.FontSize = 15;
    colormap([1 0 0;      %// red
      0 1 0;      %// green
      0 0 1;      %// blue
      .5 .5 .5])  %// grey
    hText = findobj(p,'Type','text'); % text object handles
    percentValues = get(hText,'String'); % percent values
    combinedtxt = strcat(label,percentValues'); % strings and percent values
    hText(1).String = combinedtxt(1);
    hText(2).String = combinedtxt(2);
    hText(3).String = combinedtxt(3);
    hText(4).String = combinedtxt(4);
    title(['Proteome partition, time = ', num2str(time_ss(i))], 'FontSize', 15);
    hold on
    axes(ax2);
    scatter(time_ss([1:10:4000, 4100:10:i]), energy_array([1:10:4000, 4100:10:i]), 'c')
    ylim([0 8*10^5])
    title('Energy', 'FontSize', 15);
    ylabel('a.u.', 'FontSize', 15);
    xlabel('Time (mins)', 'FontSize', 15);
    axes(ax3);
    scatter(time_ss([1:10:4000, 4100:10:i]), total_inst_gr_array([1:10:4000, 4100:10:i]), 'm')
    title('Inst. Growth Rate', 'FontSize', 15);
    ylabel('Hr^-^1', 'FontSize', 15);
    xlabel('Time (mins)', 'FontSize', 15);
    ylim([0 3])
    fig = gcf;
    fig.CurrentAxes = ax;
    %FigH = figure('Position', get(0, 'Screensize'));
    %M(idx) = getframe(FigH);
    M(idx) = getframe;
    ax.Units = 'normalized';
    disp(i)
    idx = idx+1;
    hold off
    %close
end
video = VideoWriter('04_06_ex21_2.avi', 'Uncompressed AVI');
video.FrameRate = 10;  % Default 30
open(video)
writeVideo(video, M);
close(video);

