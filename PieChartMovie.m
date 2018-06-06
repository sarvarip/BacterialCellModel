close all
load('FYP_04_06_par_ex_13')
idx = 1;
aac_array = [7500, 300, 300, 300];
for i = 1:10:4000
    %figure
    if P_count_vec_array(i, 4) == 0
        ax1 = axes('Position',[0.2 0.1 0.7 0.7],'Box','off');
        ax2 = axes('Position',[0.1 0.7 0.2 0.2],'Box','off');
        ax = axes('Position',[0 0 1 1],'Box','off');
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
        ax1 = axes('Position',[0.2 0.1 0.7 0.7],'Box','off');
        ax2 = axes('Position',[0.1 0.7 0.2 0.2],'Box','off');
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
        hText = findobj(p,'Type','text'); % text object handles
        percentValues = get(hText,'String'); % percent values
        combinedtxt = strcat(label,percentValues'); % strings and percent values
        hText(1).String = combinedtxt(1);
        hText(2).String = combinedtxt(2);
        hText(3).String = combinedtxt(3);
        t = p(8);
        t.FontSize = 15;
        colormap([1 0 0;      %// red
          0 1 0;      %// green
          0 0 1;      %// blue
          .5 .5 .5])  %// grey
        hText(4).String = combinedtxt(4);
    end
    title('Proteome partition', 'FontSize', 15);
    axes(ax2);
    bar(energy_array(i), 'c')
    ylim([0 15000])
    title('Energy', 'FontSize', 15);
    set(gca,'xticklabel','')
    fig = gcf;
    fig.CurrentAxes = ax;
    M(idx) = getframe;
    ax.Units = 'normalized';
    disp(i)
    idx = idx+1;
    hold off
    %close
end
for i = 4100:100:10000
    figure
    ax1 = axes('Position',[0.2 0.1 0.7 0.7],'Box','off');
    ax2 = axes('Position',[0.1 0.7 0.2 0.2],'Box','off');
    ax = axes('Position',[0 0 1 1],'Box','off');
    axes(ax1)
    label = {'R-proteins', 'Enzyme', 'Housekeeping', 'Heterologous'};
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
    title('Proteome partition', 'FontSize', 15);
    axes(ax2);
    bar(energy_array(i), 'c')
    ylim([0 15000])
    title('Energy', 'FontSize', 15);
    set(gca,'xticklabel','')
    fig = gcf;
    fig.CurrentAxes = ax;
    M(idx) = getframe;
    ax.Units = 'normalized';
    disp(i)
    idx = idx+1;
    hold off
    close
end
video = VideoWriter('04_06_ex13_2.avi', 'Uncompressed AVI');
video.FrameRate = 15;  % Default 30
open(video)
writeVideo(video, M);
close(video);

