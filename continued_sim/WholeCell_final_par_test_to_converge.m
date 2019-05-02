% © Peter Sarvari, University of Southern California

%issues: would be nice to save the variable time and maybe run until a
%certain time is achieved rather than up to a certain number of iterations

clc, clear

parpool(3) 
index_array = [5 7 8];
files = ["FYP_04_06_par_ex_5.mat", "FYP_04_06_par_ex_7.mat", "FYP_04_06_par_ex_8.mat"];

parfor loopvar=1:numel(index_array) 

    %Initialization from previous simulation result
    s = load(files(loopvar), 'P_count_vec_array', 'S_i', 'energy', 'location_array', 'state_array', 'temp', 'total_transcript', 'type_idx_array');     
    
    state_array = s.state_array;
    location_array = s.location_array;
    energy = s.energy;
    S_i = s.S_i;
    P_count_vec = s.P_count_vec_array(end,:);
    total_transcript = s.total_transcript;
    temp = s.temp;
    type_idx_array = s.type_idx_array;
    
    index = index_array(loopvar);
    
    no_types_mRNA = 4;
    aac_array = [7500, 300, 300, 300];
    max_elongation = 126;
    R0 = P_count_vec(1); 

    ss_start = 18000; %9000
    ss_end = 20000; %10000
    maxsteps = 10000*ss_end; 
    ref = 10000*ss_start;

    het_init_array = [0.1,1,10];
    het_rel_trc_array = [0.1,1,10];

    codon_choice = zeros(3,30);
    codon_choice(1,:) = max_elongation*ones(1,30);
    codon_choice(2,:) = max_elongation*[ones(1,9), 0.005, ones(1,20)];
    codon_choice(3,:) = max_elongation*[ones(1,20), 0.005, ones(1,9)];

    base3 = dec2base(index,3,3);
    j = str2double(base3(3))+1; %transcription rate option - function call
    k = str2double(base3(2))+1; %initiation rate option - het_init_array index
    l = str2double(base3(1))+1; %codon options - codon choice array row index

    betas = cell(1,no_types_mRNA);
    init_rates = [1,1,1,het_init_array(k)];
    betas{1} = [init_rates(1), max_elongation*ones(1,750)]; %R
    betas{2} = [init_rates(2), max_elongation*ones(1,30)]; %E
    betas{3} = [init_rates(3), max_elongation*ones(1,30)]; %Q
    betas{4} = [init_rates(4), codon_choice(l,:)]; %H 

    time = zeros(1,maxsteps+1);
    time_P_cell = cell(1,no_types_mRNA);

    %Iterations
    
    %one extra call without last argument so that we calculate
    %transition_array (cannot load it because it was not saved)
    [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(2), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_final_par(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(1), time_P_cell, P_count_vec, temp, het_rel_trc_array(j));

    exists_reference = 0;
    
    total_inst_gr_array=zeros(1, ss_end+1);
    time_ss = zeros(1, ss_end+1);
    
    P_count_vec_array=zeros(ss_end+1, 4);
    P_count_vec_array(1,:) = P_count_vec;
    energy_array = zeros(1, ss_end+1);
    energy_array(1) = energy;
    
    tic
    for timestep=2:maxsteps
        
        [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_final_par(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, het_rel_trc_array(j), transition_array);
        %reconstruct the transcripts in case changes happened
        
        if rem(timestep,10000)==0
                if exists_reference == 0
                    time_ref = time(timestep+1);
                    mass_ref = zeros(1, no_types_mRNA);
                    for i=1:no_types_mRNA-1
                        %-1, because heterologous protein should not count in mass
                        %production, since it has no benefit to the cell and it
                        %might be taken away from the cell constantly
                        mass_ref(i) = aac_array(i)*length(time_P_cell{i});
                    end
                    total_mass_ref = sum(mass_ref);
                    exists_reference = 1;
                    P_count_vec_array((timestep/10000)+1,:) = P_count_vec;
                    energy_array((timestep/10000)+1) = energy;
                else
                    time_current = time(timestep+1);
                    time_elapsed = time_current - time_ref;
                    mass_current = zeros(1, no_types_mRNA);
                    for i=1:no_types_mRNA-1
                        mass_current(i) = aac_array(i)*length(time_P_cell{i});
                    end
                    %mass_change = mass_current - mass_ref;
                    %inst_gr = (mass_change./mass_ref)/time_elapsed;
                    total_mass_current = sum(mass_current);
                    total_mass_change = total_mass_current - total_mass_ref;
                    total_inst_gr = total_mass_change/(10^8*time_elapsed)*60;
                    time_ref = time_current;
                    total_mass_ref = total_mass_current;
                    %inst_gr_array = [inst_gr_array; inst_gr]; %waste calculating both here, second can be calc from first
                    disp(['Instantaneous growth rate: ',num2str(total_inst_gr)]);
                    P_count_vec_array((timestep/10000)+1,:) = P_count_vec;
                    energy_array((timestep/10000)+1) = energy;
                    time_ss((timestep/10000)+1) = time_ref;
                    total_inst_gr_array((timestep/10000)+1) = total_inst_gr;
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
    growth_rate = (aac_array(1:3)*(total_P(1:3)-transient_P(1:3))')/(10^8*time_elapsed)*60;
    disp(['Total Growth rate: ', num2str(growth_rate)]);

    %Definition: 
    % log_no_new = log((aac_array*total_P')/10^8 +1);
    % log_no_old = log((aac_array*transient_P')/10^8 +1);
    % no_generations = (log_no_new-log_no_old)/log(2);
    % generation_time = (time_elapsed/60)/no_generations;
    % growth_rate = 1/generation_time;

    %mathematically the same as 
    % log((aac_array*total_P'+10^8)/(aac_array*transient_P'+10^8))/(log(2)*time_elapsed)*60

    %Now we only model one cell, so log_no_old = 0 (no_old is always one)
    %However, still we cannot use this, since we only model one cell and the
    %logarithmic model assumes that as soon as you have more than one, it helps
    %producing more mass/cells, but here we only model a single cell

    %approx.
    %(aac_array*total_P'-aac_array*transient_P')/(10^8*time_elapsed)*60

    avg_inst_growth_rate = mean(total_inst_gr_array(ss_start+1:end));
    std_inst_growth_rate = std(total_inst_gr_array(ss_start+1:end));
    disp(['Avg. Inst. Growth rate: ', num2str(avg_inst_growth_rate)]);
    disp(['Std. Inst. Growth rate: ', num2str(std_inst_growth_rate)]);
    
    str = ['FYP_2019_converge_par_ex_', num2str(index)];
    parsave_final(str, time_ss, state_array, location_array, type_idx_array, temp, energy, energy_array, S_i, total_transcript, total_inst_gr_array, P_count_vec_array, total_P, transient_P, P_ss, time_elapsed, production_rate, growth_rate, avg_inst_growth_rate, std_inst_growth_rate);
    
end