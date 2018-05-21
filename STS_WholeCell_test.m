% © Peter Sarvari, Imperial College London

%Initialization
total_mRNA = 3971;
no_types_mRNA = 3;

total_protein_init = [7373, 10924, 138080];
total_transcript = [187, 278, 3506];
temp = [1*ones(1,total_transcript(1)), 1*ones(1,total_transcript(2)), 1*ones(1,total_transcript(3))]; 
init_rates = [1,1,1];
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

energy = 8.05*10^6;
S_i = 128;
R0 = total_protein_init(1); 

maxsteps = 2000*36000; %~est. 10 hrs
ref = 2000*25000;
time = zeros(1,maxsteps+1);
time_P_cell = []; %CHANGED! ONLY A VECTOR, since whenever a protein is produced, so mass changes, all protein amounts change a bit due to dilution
P_count_vec = total_protein_init;
state_array = ones(1,R0);
location_array = zeros(1,R0);

%Iterations

tic
for timestep=1:maxsteps
    if timestep == 1
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_WholeCell(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp);
        %reconstruct the transcripts in case changes happened
    elseif timestep == ref+1 
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_WholeCell(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
        tot_mass_ref = [7500, 300, 300]*P_count_vec';
    elseif timestep == maxsteps
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_WholeCell(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
        tot_mass_end = [7500, 300, 300]*P_count_vec';
    else    
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_WholeCell(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
        %reconstruct the transcripts in case changes happened
    end
    
    if rem(timestep,2000)==0
        ratio = timestep/maxsteps;
        disp([num2str(100*ratio), '%']);
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
disp(production_rate);
disp(growth_rate);