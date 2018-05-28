% © Peter Sarvari, Imperial College London

clc, clear

%parpool(4) %test - works, tested!
parpool(14) %run1
%parpool(13) %run2

parfor index=1:14 %run 1 %1:4 for test
%parfor index=15:27 %run2

    %Heterologous version
    %Initialization
    total_mRNA = 3971;
    no_types_mRNA = 4;
    total_protein_init = [7373, 10924, 138080, 0];
    total_transcript = [187, 278, 3506, 0];
    aac_array = [7500, 300, 300, 300];
    max_elongation = 126;
    energy = 10^5;
    S_i = 128;
    R0 = total_protein_init(1); 

    het_init_array = [0.1,1,10];

    codon_choice = zeros(3,30);
    codon_choice(1,:) = max_elongation*ones(1,30);
    codon_choice(2,:) = max_elongation*[ones(1,9), 0.1, ones(1,20)];
    codon_choice(3,:) = max_elongation*[ones(1,20), 0.1, ones(1,9)];

    type_idx_array = [ones(1,total_transcript(1)), 2*ones(1,total_transcript(2)), 3*ones(1,total_transcript(3)), 4*ones(1,total_transcript(4))]; 

    ss_start = 25000;
    ss_end = 36000;
    maxsteps = 2000*ss_end; %~est. 10 hrs
    ref = 2000*ss_start;

    base3 = dec2base(index,3,3);
    j = str2double(base3(1))+1; %transcription rate option - function call
    k = str2double(base3(2))+1; %initiation rate option - het_init_array index
    l = str2double(base3(3))+1; %codon options - codon choice array row index


    betas = cell(1,no_types_mRNA);
    init_rates = [1,1,1,het_init_array(k)];
    betas{1} = [init_rates(1), max_elongation*ones(1,750)]; %R
    betas{2} = [init_rates(2), max_elongation*ones(1,30)]; %E
    betas{3} = [init_rates(3), max_elongation*ones(1,30)]; %Q
    betas{4} = [init_rates(4), codon_choice(l,:)]; %H % INPUT! Change accordingly!
    temp = [init_rates(1)*ones(1,total_transcript(1)), init_rates(2)*ones(1,total_transcript(2)), init_rates(3)*ones(1,total_transcript(3)), init_rates(4)*ones(1,total_transcript(4))]; 

    time = zeros(1,maxsteps+1);
    time_P_cell = cell(1,no_types_mRNA);
    P_count_vec = total_protein_init;
    state_array = ones(1,R0);
    location_array = zeros(1,R0);

    %Iterations

    exists_reference = 0;
    total_inst_gr_array=zeros(1, ss_end-ss_start);
    P_count_vec_array=zeros(ss_end-ss_start+1, 4);
    time_ss = zeros(1, ss_end-ss_start+1);
    tic
    
    if j == 1
    
        for timestep=1:maxsteps
            if timestep == 1
                [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final_ex1(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp);
                %reconstruct the transcripts in case changes happened
            else    
                [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final_ex1(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
                %reconstruct the transcripts in case changes happened
            end

            if rem(timestep,2000)==0
                if timestep >= ref
                    if exists_reference == 0
                        time_ref = time(timestep+1);
                        mass_ref = zeros(1, no_types_mRNA); 
                        for i=1:no_types_mRNA
                            mass_ref(i) = aac_array(i)*length(time_P_cell{i});
                        end
                        total_mass_ref = sum(mass_ref);
                        exists_reference = 1;
                        P_count_vec_array((timestep/2000)-ss_start+1,:) = P_count_vec;
                        time_ss((timestep/2000)-ss_start+1) = time_ref;
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
                        total_inst_gr = total_mass_change/(total_mass_ref*time_elapsed);
                        time_ref = time_current;
                        mass_ref = mass_current;
                        %inst_gr_array = [inst_gr_array; inst_gr]; %waste calculating both here, second can be calc from first
                        total_inst_gr_array((timestep/2000)-ss_start) = total_inst_gr;
                        disp(['Instantaneous growth rate: ',num2str(inst_gr)]);
                        P_count_vec_array((timestep/2000)-ss_start+1,:) = P_count_vec;
                        time_ss((timestep/2000)-ss_start+1) = time_ref;
                    end
                end    
                ratio = timestep/maxsteps;
                disp([num2str(100*ratio), '%']);
                disp(['Energy: ',num2str(energy)]);
            end
        end
        toc


    elseif j == 2
        
        for timestep=1:maxsteps
            if timestep == 1
                [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final_ex2(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp);
                %reconstruct the transcripts in case changes happened
            else    
                [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final_ex2(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
                %reconstruct the transcripts in case changes happened
            end

            if rem(timestep,2000)==0
                if timestep >= ref
                    if exists_reference == 0
                        time_ref = time(timestep+1);
                        mass_ref = zeros(1, no_types_mRNA); 
                        for i=1:no_types_mRNA
                            mass_ref(i) = aac_array(i)*length(time_P_cell{i});
                        end
                        total_mass_ref = sum(mass_ref);
                        exists_reference = 1;
                        P_count_vec_array((timestep/2000)-ss_start+1,:) = P_count_vec;
                        time_ss((timestep/2000)-ss_start+1) = time_ref;
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
                        total_inst_gr = total_mass_change/(total_mass_ref*time_elapsed);
                        time_ref = time_current;
                        mass_ref = mass_current;
                        %inst_gr_array = [inst_gr_array; inst_gr]; %waste calculating both here, second can be calc from first
                        total_inst_gr_array((timestep/2000)-ss_start) = total_inst_gr;
                        disp(['Instantaneous growth rate: ',num2str(inst_gr)]);
                        P_count_vec_array((timestep/2000)-ss_start+1,:) = P_count_vec;
                        time_ss((timestep/2000)-ss_start+1) = time_ref;
                    end
                end    
                ratio = timestep/maxsteps;
                disp([num2str(100*ratio), '%']);
                disp(['Energy: ',num2str(energy)]);
            end
        end
        toc
        
    else %j==3
        
        for timestep=1:maxsteps
            if timestep == 1
                [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final_ex3(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp);
                %reconstruct the transcripts in case changes happened
            else    
                [state_array, location_array, type_idx_array, total_transcript, energy, S_i, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = WholeCell_2_final_ex3(state_array, location_array, betas, type_idx_array, total_transcript, energy, S_i, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
                %reconstruct the transcripts in case changes happened
            end

            if rem(timestep,2000)==0
                if timestep >= ref
                    if exists_reference == 0
                        time_ref = time(timestep+1);
                        mass_ref = zeros(1, no_types_mRNA); 
                        for i=1:no_types_mRNA
                            mass_ref(i) = aac_array(i)*length(time_P_cell{i});
                        end
                        total_mass_ref = sum(mass_ref);
                        exists_reference = 1;
                        P_count_vec_array((timestep/2000)-ss_start+1,:) = P_count_vec;
                        time_ss((timestep/2000)-ss_start+1) = time_ref;
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
                        total_inst_gr = total_mass_change/(total_mass_ref*time_elapsed);
                        time_ref = time_current;
                        mass_ref = mass_current;
                        %inst_gr_array = [inst_gr_array; inst_gr]; %waste calculating both here, second can be calc from first
                        total_inst_gr_array((timestep/2000)-ss_start) = total_inst_gr;
                        disp(['Instantaneous growth rate: ',num2str(inst_gr)]);
                        P_count_vec_array((timestep/2000)-ss_start+1,:) = P_count_vec;
                        time_ss((timestep/2000)-ss_start+1) = time_ref;
                    end
                end    
                ratio = timestep/maxsteps;
                disp([num2str(100*ratio), '%']);
                disp(['Energy: ',num2str(energy)]);
            end
        end
        toc
        
    end
    
    %Final average production calcs

    time_elapsed = time(end)-time(ref);
    transient_P=zeros(1,no_types_mRNA);
    total_P = zeros(1,no_types_mRNA);
    for i=1:no_types_mRNA
        transient_P(i) = sum(time_P_cell{i}<=time(ref)) + total_protein_init(i); %we only care about states from the 1001th timstep on (incl), but
        total_P(i) = length(time_P_cell{i});
        %if a protein was produced at time = 1001, it was the trasition from state
        %1000 to 1001, so we disregard it! CHANGED: we are not averaging
        %here; also ss_start is estimated, for simplicity, use ref instead
        %of ref+1!!
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
    
    
    str = ['FYP_29_05_par_ex_', num2str(index)];
    parsave(str, total_inst_gr_array, P_count_vec_array, time_ss, total_P, transient_P, P_ss, time_elapsed, production_rate, growth_rate, avg_inst_growth_rate, std_inst_growth_rate);
    
end