% © Peter Sarvari, Imperial College London

clc, clear

%Heterologous version
%Initialization
total_mRNA = 3971;
no_types_mRNA = 4;

total_protein_init = [7373, 10924, 138080, 0];
total_transcript = [187, 278, 3506, 0];
init_rates = [1,1,1,0.1];
temp = [init_rates(1)*ones(1,total_transcript(1)), init_rates(2)*ones(1,total_transcript(2)), init_rates(3)*ones(1,total_transcript(3)), init_rates(4)*ones(1,total_transcript(4))]; 
max_elongation = 126;
het_mass = 300; % INPUT! Change accordingly!

% cum_total_transcript = cumsum(total_transcript);
% betas = cell(1,total_mRNA); %TO DO: SHOULD BE JUST no_types_mRNA here as well!!
% 
% type1_rna_mat = [1, 126*ones(1,750)];
% betas(1:total_transcript(1)) = {type1_rna_mat}; %R, 187 copy
% 
% type2_rna_mat = [1, 126*ones(1,30)];
% betas((total_transcript(1)+1):cum_total_transcript(2)) = {type2_rna_mat}; %E, 278 copy
% 
% type3_rna_mat = [1, 126*ones(1,30)];
% betas((total_transcript(2)+1):cum_total_transcript(3)) = {type3_rna_mat}; %Q, 3506 copy

betas = cell(1,no_types_mRNA);
betas{1} = [init_rates(1), max_elongation*ones(1,750)]; %R
betas{2} = [init_rates(2), max_elongation*ones(1,30)]; %E
betas{3} = [init_rates(3), max_elongation*ones(1,30)]; %Q
betas{4} = [init_rates(4), max_elongation*ones(1,het_mass/10)]; %H % INPUT! Change accordingly!
type_idx_array = [ones(1,total_transcript(1)), 2*ones(1,total_transcript(2)), 3*ones(1,total_transcript(3)), 4*ones(1,total_transcript(4))]; 

energy = 10^5;
S_i = 128;
R0 = total_protein_init(1); 

maxsteps = 2000*36000; %~est. 10 hrs
ref = 2000*25000;
time = zeros(1,maxsteps+1);
time_P_cell = cell(1,no_types_mRNA);
P_count_vec = total_protein_init;
state_array = ones(1,R0);
location_array = zeros(1,R0);

%Iterations

exists_reference = 0;
inst_gr_array = [];
tic
for timestep=1:maxsteps
    if timestep == 1
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final_ex1(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp);
        %reconstruct the transcripts in case changes happened
    elseif timestep == ref+1 
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final__ex1(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
        tot_mass_ref = [7500, 300, 300, het_mass]*P_count_vec';
    elseif timestep == maxsteps
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final__ex1(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
        tot_mass_end = [7500, 300, 300, het_mass]*P_count_vec';
    else    
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final__ex1(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
        %reconstruct the transcripts in case changes happened
    end
    
    if rem(timestep,2000)==0
        if timestep >= ref+1
            if exists_reference == 0
                time_ref = time(timestep+1);
                mass_ref = [7500, 300, 300, het_mass]*P_count_vec';
                exists_reference = 1;
            else
                time_current = time(timestep+1);
                time_elapsed = time_current - time_ref;
                mass_current = [7500, 300, 300, het_mass]*P_count_vec';
                mass_change = mass_current - mass_ref;
                inst_gr = mass_change/time_elapsed;
                time_ref = time_current;
                mass_ref = mass_current;
                inst_gr_array = [inst_gr_array, inst_gr];
                disp(['Instantaneous growth rate: ',num2str(inst_gr)]);
            end
        end    
        ratio = timestep/maxsteps;
        disp([num2str(100*ratio), '%']);
        disp(['Energy: ',num2str(energy)]);
    end
end
toc

%Final average production calcs

time_elapsed = time(end)-time(ref+1);
transient_P=zeros(1,no_types_mRNA);
for i=1:no_types_mRNA
    transient_P(i) = sum(time_P_cell{i}<=time(ref+1)) + total_protein_init; %we only can about states from the 1001th timstep on (incl), but
    %if a protein was produced at time = 1001, it was the trasition from state
    %1000 to 1001, so we disregard it!
end

P_ss = P_count_vec-transient_P;
production_rate = P_ss/time_elapsed;
growth_rate = (tot_mass_end-tot_mass_ref)/time_elapsed;
disp(['Production rate: ', num2str(production_rate)]);
disp(['Total Growth rate: ', num2str(growth_rate)]);
avg_inst_growth_rate = mean(inst_gr_array);
std_inst_growth_rate = std(inst_gr_array);
disp(['Avg. Inst. Growth rate: ', num2str(avg_inst_growth_rate)]);
disp(['Std. Inst. Growth rate: ', num2str(std_inst_growth_rate)]);

save('FYP_ex1_25_05')