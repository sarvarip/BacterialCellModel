close all
load('FYP_04_06_par_ex_13')
idx = 1;
aac_array = [7500, 300, 300, 300];
for i = 1:10:4000
    %figure
    ax1 = axes('Position',[0.2 0.1 0.7 0.7],'Box','off');
    ax2 = axes('Position',[0.1 0.7 0.2 0.2],'Box','off');
    ax = axes('Position',[0 0 1 1],'Box','off');
    axes(ax1)
    label = {'R-proteins', 'Enzyme', 'Housekeeping', 'Heterologous'};
    p= pie(P_count_vec_array(i, :).*aac_array, label);
    t = p(2);
    t.FontSize = 15;
    t = p(4);
    t.FontSize = 15;
    t = p(6);
    t.FontSize = 15;
    colormap([1 0 0;      %// red
          0 1 0;      %// green
          0 0 1]);      %// blue
    if P_count_vec_array(i, :) > 0
        t = p(8);
        t.FontSize = 15;
        colormap([1 0 0;      %// red
          0 1 0;      %// green
          0 0 1;      %// blue
          .5 .5 .5])  %// grey
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
    p= pie(P_count_vec_array(i, :).*aac_array, label);
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
    if P_count_vec_array(i, :) > 0
        t = p(8);
        t.FontSize = 15;
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
    close
end
video = VideoWriter('04_06_ex13.avi', 'Uncompressed AVI');
open(video)
writeVideo(video, M);
close(video);

