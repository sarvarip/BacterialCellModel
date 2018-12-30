% © Peter Sarvari, Imperial College London

clc, clear

%Initialization

load('FYP_7_06_endo_gr_newss_ns_600', 'P_count_vec', 'S_i', 'aac_array', 'betas', 'energy', 'init_rates', 'location_array', 'max_elongation', 'no_types_mRNA', 'state_array', 'temp', 'total_mRNA', 'total_transcript', 'transition_array', 'type_idx_array'); 
ss_start = 1; %45000
ss_end = 10000; %50000
maxsteps = 2000*ss_end; %~est. 20 hrs
ref = 2000*ss_start;
n_s = 600;

%Iterations

exists_reference = 0;
total_inst_gr_array=zeros(1, ss_end-ss_start+1);
P_count_vec_array=zeros(ss_end-ss_start+1, 3);
transcript_array=zeros(ss_end-ss_start+1, 3);
time_ss = zeros(1, ss_end-ss_start+1);
time = zeros(1,maxsteps+1);
time_P_cell = cell(1,no_types_mRNA);
%load('FYP_30_05_endo_newss.mat', 'time_P_cell')
%If cotinued, copy time_P_cell for growth rate purposes

tic
for timestep=1:maxsteps   
    [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_WholeCell_2_final_ns_10(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array, n_s);
    if rem(timestep,2000)==0
        if exists_reference == 0
            time_ref = time(timestep+1);
            mass_ref = zeros(1, no_types_mRNA);
            for i=1:no_types_mRNA
                mass_ref(i) = aac_array(i)*length(time_P_cell{i});
            end
            total_mass_ref = sum(mass_ref);
            exists_reference = 1;
            %time_ss((timestep/2000)-ss_start+1) = time_ref;
        else
            time_current = time(timestep+1);
            time_elapsed = time_current - time_ref;
            mass_current = zeros(1, no_types_mRNA);
            for i=1:no_types_mRNA
                mass_current(i) = aac_array(i)*length(time_P_cell{i});
            end
            %mass_change = mass_current - mass_ref;
            %inst_gr = (mass_change./mass_ref)/time_elapsed;
            total_mass_current = sum(mass_current);
            total_mass_change = total_mass_current - total_mass_ref;
            %Always divide by 10^8, since one cell is causing the
            %change to happen!
            total_inst_gr = total_mass_change/(10^8*time_elapsed)*60;
            time_ref = time_current;
            total_mass_ref = total_mass_current;
            %inst_gr_array = [inst_gr_array; inst_gr]; %waste calculating both here, second can be calc from first
            disp(['Instantaneous growth rate: ',num2str(total_inst_gr)]);
            if timestep >= ref
                time_ss((timestep/2000)-ss_start+1) = time_ref;
                total_inst_gr_array((timestep/2000)-ss_start+1) = total_inst_gr;
                P_count_vec_array((timestep/2000)-ss_start+1,:) = P_count_vec;
                transcript_array((timestep/2000)-ss_start+1,:) = total_transcript;
            end
        end
        ratio = timestep/maxsteps;
        disp([num2str(100*ratio), '%']);
        disp(['Energy: ',num2str(energy)]);
    end
end
toc

%Final average production calcs

time_elapsed = time(end)-time(ref);
transient_P=zeros(1,no_types_mRNA);
total_P = zeros(1,no_types_mRNA);
for i=1:no_types_mRNA
    transient_P(i) = sum(time_P_cell{i}<=time(ref)); %we only care about states from the 1001th timstep on (incl), but
    total_P(i) = length(time_P_cell{i});
    %if a protein was produced at time = 1001, it was the trasition from state
    %1000 to 1001, so we disregard it! CHANGED: we are not averaging
    %here; also ss_start is estimated, for simplicity, use ref instead
    %of ref+1!!
end


P_ss = total_P-transient_P;
production_rate = P_ss/time_elapsed;
disp(['Production rate: ', num2str(production_rate)]);
growth_rate = (aac_array*total_P'-aac_array*transient_P')/(10^8*time_elapsed)*60;
disp(['Total Growth rate: ', num2str(growth_rate)]);

avg_inst_growth_rate = mean(total_inst_gr_array);
std_inst_growth_rate = std(total_inst_gr_array);
disp(['Avg. Inst. Growth rate: ', num2str(avg_inst_growth_rate)]);
disp(['Std. Inst. Growth rate: ', num2str(std_inst_growth_rate)]);

save('FYP_22_12_endo_gr_newss_ns_600_transcript_arr')


%Piechart

%clear
%load('FYP_1_06_endo_gr_newss_ns_10')
%change title!

% close all
% figure;
% ax1 = axes('Position',[0.2 0.1 0.7 0.7],'Box','off');
% ax2 = axes('Position',[0.1 0.7 0.2 0.2],'Box','off');
% ax3 = axes('Position',[0.1 0.1 0.2 0.2],'Box','off');
% ax = axes('Position',[0 0 1 1],'Box','off');
% axes(ax1)
% label = {'R-proteins: ', 'Enzyme: ', 'Housekeeping: '};
% p= pie(P_count_vec_array(end, :).*aac_array);
% t = p(2);
% t.FontSize = 15;
% t = p(4);
% t.FontSize = 15;
% t = p(6);
% t.FontSize = 15;
% colormap([1 0 0;      %// red
%       0 1 0;      %// green
%       0 0 1]);      %// blue
% hText = findobj(p,'Type','text'); % text object handles
% percentValues = get(hText,'String'); % percent values
% combinedtxt = strcat(label,percentValues'); % strings and percent values
% hText(1).String = combinedtxt(1);
% hText(2).String = combinedtxt(2);
% hText(3).String = combinedtxt(3);
% title(['Endogenous Proteome partition n_s=',num2str(n_s)] , 'FontSize', 15); %CHANGE n_s!
% axes(ax2);
% bar((energy), 'c')
% title('Energy', 'FontSize', 15);
% ylabel('Arb. unit', 'FontSize', 15);
% set(gca,'xticklabel','')
% axes(ax3);
% bar(([avg_inst_growth_rate, std_inst_growth_rate]), 'm')
% title('Inst. Growth Rate', 'FontSize', 15);
% set(gca,'xticklabel',{'Average', 'Std'}, 'FontSize', 15)
% ylabel('h^-^1', 'FontSize', 15);
% fig = gcf;
% fig.CurrentAxes = ax;
% ax.Units = 'normalized';
% hold off