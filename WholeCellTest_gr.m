% © Peter Sarvari, Imperial College London

%Initialization
total_mRNA = 3971;
no_types_mRNA = 3;

total_protein_init = [7373, 10924, 138080];
total_transcript = [187, 278, 3506];
init_rates = [1,1,1];
temp = [init_rates(1)*ones(1,total_transcript(1)), init_rates(2)*ones(1,total_transcript(2)), init_rates(3)*ones(1,total_transcript(3))]; 
max_elongation = 126;

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
type_idx_array = [ones(1,total_transcript(1)), 2*ones(1,total_transcript(2)), 3*ones(1,total_transcript(3))]; 

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
total_inst_gr_array=[];
aac_array = [7500, 300, 300];
tic
for timestep=1:maxsteps
    if timestep == 1
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_WholeCell_2_final(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp);
        %reconstruct the transcripts in case changes happened
    else    
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_WholeCell_2_final(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
        %reconstruct the transcripts in case changes happened
    end
    
    if rem(timestep,2000)==0
        if timestep >= ref+1
            if exists_reference == 0
                time_ref = time(timestep+1);
                mass_ref = zeros(1, no_types_mrna); 
                for i=1:no_types_mRNA
                    mass_ref(i) = aac_array(i)*length(time_P_cell{i});
                end
                total_mass_ref = sum(mass_ref);
                exists_reference = 1;
            else
                time_current = time(timestep+1);
                time_elapsed = time_current - time_ref;
                mass_current = zeros(1, no_types_mrna); 
                for i=1:no_types_mRNA
                    mass_current(i) = aac_array(i)*length(time_P_cell{i});
                end
                %mass_change = mass_current - mass_ref;
                %inst_gr = (mass_change./mass_ref)/time_elapsed;
                total_mass_current = sum(mass_current);
                total_mass_change = total_mass_current - total_mass_ref;
                total_inst_gr = total_mass_change/(total_mass_ref*time_elapsed);
                time_ref = time_current;
                mass_ref = mass_current;
                %inst_gr_array = [inst_gr_array; inst_gr]; %waste calculating both here, second can be calc from first
                total_inst_gr_array = [total_inst_gr_array, total_inst_gr];
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
total_P = zeros(1,no_types_mRNA);
for i=1:no_types_mRNA
    transient_P(i) = sum(time_P_cell{i}<=time(ref+1)) + total_protein_init(i); %we only care about states from the 1001th timstep on (incl), but
    total_P(i) = length(time_P_cell{i});
    %if a protein was produced at time = 1001, it was the trasition from state
    %1000 to 1001, so we disregard it!
end

P_ss = total_P-transient_P;
production_rate = P_ss/time_elapsed;
disp(['Production rate: ', num2str(production_rate)]);
growth_rate = (aac_array*total_P')/((aac_array*transient_P')*time_elapsed);
disp(['Total Growth rate: ', num2str(growth_rate)]);
avg_inst_growth_rate = mean(total_inst_gr_array);
std_inst_growth_rate = std(total_inst_gr_array);
disp(['Avg. Inst. Growth rate: ', num2str(avg_inst_growth_rate)]);
disp(['Std. Inst. Growth rate: ', num2str(std_inst_growth_rate)]);

save('FYP_26_05_endo')