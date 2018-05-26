function [state_array, location_array, type_idx_array, frac_total_transcript, energy, s_i, time_out, time_P_cell, frac_P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_WholeCell_2_final (state_array, location_array, beta_array_cell, type_idx_array, frac_total_transcript, energy, s_i, time_in, time_P_cell, frac_P_count_vec, temp, transition_array)
    
    % � Peter Sarvari, Imperial College London

    n_s = 100; % INPUT! Change accordingly!
    
    total_transcript = round(frac_total_transcript);
    P_count_vec = round(frac_P_count_vec);
    
    R = length(state_array);
    no_mrna = length(type_idx_array); 
    ribo_movement = 0;
    protein_diluted = 0;

    if nargin == 11
        transition_array = zeros(no_mrna,R);
        for mrna = 1:no_mrna
            transition_array(mrna,:) = beta_array_cell{type_idx_array(mrna)}(state_array.*(location_array==mrna)+1*(location_array~=mrna)); %otherwise error because e.g. site 20 may not exist on a short mRNA, add one for places that will be zeroed in the next line, otherwise zero indexing error
            transition_array(mrna,:) = transition_array(mrna,:).*(location_array==mrna); %dot multiply so that transition array stays R long for each mRNA
        end
        transition_array = sum(transition_array, 1);
    end
    
    free = sum(location_array==0);
    
    %R transcription rate
    R_tr = 27*energy/(energy+30000); %430
    %E transcription rate
    E_tr = 20*energy/(energy+3000); %4.4
    %Q transcription rate
    Q_tr = 850*(energy/(energy+10000))*((5000^4)/(5000^4+frac_total_transcript(3)^4)); %4.4, 152000
    %RNA degradation with ribosome protection
    unoccupied_mrna = setdiff(1:no_mrna, unique(location_array));
    no_unoccupied_mrna = length(unoccupied_mrna);
    if no_unoccupied_mrna > 0
        deg = 0.1; %might wanna change to 0.2 to ensure smaller number of transcripts 
        %and compensate for ribosome protection; moreover bionumbers also suggests that
    else 
        deg = 0;
    end
    %Nutrient transport
    transport = 0.363*frac_P_count_vec(2);
    %Nutrient metabolism
    met = 2.9*frac_P_count_vec(2)*s_i^3/(100^3+s_i^3);
    
    if energy >= 1
        %elongation rate
        %need to make sure that conditions are such that at energy=1 the
        %elongation rates are small; otherwise back and forth oscillations!
        translation_scale = (energy^2/(energy^2+10000^2));
        transition_array(transition_array>0) = transition_array(transition_array>0)*translation_scale; 
        %division at the end to go back to normal so that due to recursion the
        %multiplication does not accumulate!
    
        beta_array_cell{1}(2:end) = beta_array_cell{1}(2:end)*translation_scale;
        beta_array_cell{2}(2:end) = beta_array_cell{2}(2:end)*translation_scale;
        beta_array_cell{3}(2:end) = beta_array_cell{3}(2:end)*translation_scale;
        
        transition_array = [sum(temp)*free, transition_array, R_tr, E_tr, Q_tr, deg*sum(frac_total_transcript), transport, met];
    else
        transition_array_to_return = transition_array;
        temporary_transition_array = zeros(1,R);
        transition_array = [sum(temp)*free, temporary_transition_array, R_tr, E_tr, Q_tr, deg*sum(frac_total_transcript), transport, met];
    end
    
    
    total_out_rate = sum(transition_array);
    
    deltaT = exprnd(1/total_out_rate);
    time_out = time_in + deltaT;
    
    saved_rand = rand;
    
    idx = find(cumsum(transition_array/total_out_rate)>saved_rand, 1, 'first');
    tt = 1:R;
    free_ribos = tt(location_array==0);
    if idx == length(transition_array) %metabolism
        s_i = s_i-1;
        energy = energy+n_s;
        %disp('metabolism');
    elseif idx == length(transition_array)-1 %transport
        s_i = s_i+1;
        %disp('Transport');
    elseif idx == length(transition_array)-2 %mRNA degradation
        unoccupied_types = unique(type_idx_array(unoccupied_mrna));
        tt = frac_total_transcript(unoccupied_types);
        saved_rand3 = rand;
        idx3 = find(cumsum(tt/sum(tt))>saved_rand3, 1, 'first');
        possible_degrade = unoccupied_mrna(type_idx_array(unoccupied_mrna)==unoccupied_types(idx3));
        to_degrade = possible_degrade(randi(length(possible_degrade)));
        %deletion from database
        total_transcript(unoccupied_types(idx3)) = total_transcript(unoccupied_types(idx3))-1;
        frac_total_transcript(unoccupied_types(idx3)) = frac_total_transcript(unoccupied_types(idx3))-1;
        type_idx_array(to_degrade) = []; 
        location_array = location_array - (location_array>to_degrade); %shift everything back after the 
        %degraded mRNA to that numbering is continuous; otherwise
        %type_idx_array size decreases but maximum in location_array is
        %still the same! 
        no_mrna = no_mrna - 1;
        temp(to_degrade) = []; 
        %disp('mRNA degraded');
        if max(location_array) > no_mrna
            disp('Indexing Error with location_array! Program paused, Press Ctrl-C to exit!');
            pause;
        end
    elseif idx == length(transition_array)-3 %Q type mRNA trascription
        total_transcript(3) = total_transcript(3)+1;
        frac_total_transcript(3) = frac_total_transcript(3)+1;
        type_idx_array = [type_idx_array, 3];
        temp = [temp, beta_array_cell{3}(1)]; 
        %disp('Q mRNA produced');
    elseif idx == length(transition_array)-4 %E type mRNA trascription
        total_transcript(2) = total_transcript(2)+1;
        frac_total_transcript(2) = frac_total_transcript(2)+1;
        type_idx_array = [type_idx_array, 2];
        temp = [temp, beta_array_cell{2}(1)]; 
        %disp('E mRNA produced');
    elseif idx == length(transition_array)-5 %R type mRNA trascription
        total_transcript(1) = total_transcript(1)+1;
        frac_total_transcript(1) = frac_total_transcript(1)+1;
        type_idx_array = [type_idx_array, 1];
        temp = [temp, beta_array_cell{1}(1)]; 
        %disp('R mRNA produced');
    elseif idx == 1 %initiation
        saved_rand2 = rand;
        idx2 = find(cumsum(temp/sum(temp))>saved_rand2, 1, 'first'); %note that for this, 
        %temp doesn't have to be ordered, such that type 2 follows type 1, etc
        ribo_loc = idx2;
        ribo = free_ribos(randi(length(free_ribos))); 
        ribo_movement = 1;
        %disp('Initiation');
    else %elongation
        ribo = idx-1;
        ribo_loc = location_array(ribo);
        energy = energy-1;
        ribo_movement = 1;
        %disp('Elongation');
    end
    
    if ribo_movement == 1
        
        prev_state = state_array(ribo);
        prev_location = location_array(ribo);
        
        %update
        if location_array(ribo)~=0
            if state_array(ribo) == length(beta_array_cell{type_idx_array(location_array(ribo))})
                state_array(ribo) = 1;
                location_array(ribo) = 0; %transition array is reset at end
                old_mass = [7500, 300, 300]*frac_P_count_vec';
                frac_P_count_vec(type_idx_array(ribo_loc)) = frac_P_count_vec(type_idx_array(ribo_loc))+1;
                P_count_vec(type_idx_array(ribo_loc)) = P_count_vec(type_idx_array(ribo_loc))+1;
                time_P_cell{type_idx_array(ribo_loc)} = [time_P_cell{type_idx_array(ribo_loc)}, time_out];
                if type_idx_array(ribo_loc) == 1 %if R protein is produced, extend rbosome data arrays
                    state_array = [state_array, 1];
                    location_array = [location_array, 0];
                    transition_array = [transition_array(1:end-6), 0, transition_array(end-5:end)];
                    %disp('R protein produced');
                end
                new_mass = [7500, 300, 300]*frac_P_count_vec';
                ratio = new_mass/old_mass;
                %disp('Protein Produced');
                
                %dilution
                energy = energy/ratio;
                frac_total_transcript = frac_total_transcript/ratio;
                frac_P_count_vec = frac_P_count_vec/ratio;
                s_i = s_i/ratio;
                %rounding issues
                if ~isequal(round(frac_total_transcript), total_transcript)
                    decreased = round(frac_total_transcript)~=total_transcript; %logical array
                    types = 1:3;
                    decreased_types = types(decreased);
                    %deleting transcripts, which have smallest amount of ribos
                    %on them
                    %Consider: are there transcripts w/o ribosomes?
                    unoccupied_types = unique(type_idx_array(unoccupied_mrna));
                    for i = decreased_types
                        %If there are none
                        if ~ismember(i, unoccupied_types)
                            location_array_nonzero = location_array(location_array>0);
                            mrna_to_be_deleted = mode_least(location_array_nonzero(type_idx_array(location_array_nonzero)==i)); %if I leave zero elements in then indexing error
                            %reset ribosomes on that mRNA
                            ribos_on_mrna_to_be_deleted = find(location_array==mrna_to_be_deleted);
                            location_array(ribos_on_mrna_to_be_deleted) = 0;
                            state_array(ribos_on_mrna_to_be_deleted) = 1;
                            transition_array(ribos_on_mrna_to_be_deleted+1) = 0; 
                            %delete mRNA
                            type_idx_array(mrna_to_be_deleted) = [];
                            no_mrna = no_mrna - 1;
                            location_array = location_array - (location_array>mrna_to_be_deleted);
                            temp(mrna_to_be_deleted) = [];
                            %disp('non-empty mRNA diluted');
                        else
                            unoccupied_mrna = setdiff(1:no_mrna, unique(location_array)); %need to calculate again, because no_mrna 
                            %might have been reduced in previous loop!!
                            unoccupied_correct_type = unoccupied_mrna(type_idx_array(unoccupied_mrna)==i);
                            mrna_to_be_deleted = unoccupied_correct_type(randi(length(unoccupied_correct_type))); %just take first one
                            %delete mRNA
                            type_idx_array(mrna_to_be_deleted) = [];
                            no_mrna = no_mrna - 1;
                            temp(mrna_to_be_deleted) = [];
                            location_array = location_array - (location_array>mrna_to_be_deleted);
                            %disp('empty mRNA diluted');
                        end
                    end
                end
                if round(frac_P_count_vec(1)) ~= P_count_vec(1)
                    %if there are free ribosomes - there must be since one just
                    %got released, so might as well just delete that!!!
                    state_array(ribo) = [];
                    location_array(ribo) = [];
                    transition_array(ribo+1) = [];
                    %disp('Protein diluted');
                    protein_diluted = 1;
                end
            else
                state_array(ribo) = state_array(ribo) + 1;
            end
        else
            state_array(ribo) = 2;
            location_array(ribo) = ribo_loc;
        end
        
        %taking care of previous ribos
        if prev_state == 2
            temp(prev_location) = beta_array_cell{type_idx_array(prev_location)}(1);
        elseif prev_state == 1
            temp(ribo_loc) = 0;
        elseif sum(state_array==(prev_state-1))>0
            if sum(location_array(state_array==(prev_state-1))==prev_location) == 1
                vec = find(state_array==(prev_state-1));
                rib = vec(location_array(state_array==(prev_state-1))==prev_location);
                transition_array(rib+1) = beta_array_cell{type_idx_array(prev_location)}(prev_state-1);
            end
        end
        
        if protein_diluted == 0
            %taking care of updated transition rate for chosen ribo 
            %only if protein wasn't deleted (diluted)
            next_state = state_array(ribo);
            next_location = location_array(ribo);
            
            if next_location ~= 0
                if sum(state_array==(next_state+1))>0
                    if sum(location_array(state_array==(next_state+1))==next_location) == 0
                        transition_array(ribo+1) = beta_array_cell{type_idx_array(next_location)}(next_state);
                    else
                        transition_array(ribo+1) = 0;
                    end
                else
                    transition_array(ribo+1) = beta_array_cell{type_idx_array(next_location)}(next_state);
                end
            else
                transition_array(ribo+1) = 0;
            end
        end
    end
    
    %Debug I
    
    free = sum(location_array==0);
    tt = 1:length(state_array);
    free_ribos = tt(location_array==0);
    
    if sum(transition_array(free_ribos+1)==0) ~= free
        disp('Location-State-Transition Error. Program paused, Press Ctrl-C to exit!');
        pause;
    end
    
    %Debug I ends
    
    if energy==1
        if exist('transition_array_to_return', 'var')
            transition_array = transition_array_to_return;
        else
            transition_array = transition_array(2:end-6);
            transition_array(transition_array>0) = transition_array(transition_array>0)/translation_scale; 
        end
    else
        transition_array = transition_array(2:end-6);
        transition_array(transition_array>0) = transition_array(transition_array>0)/translation_scale; 
    end
    
    %Debug II
    
    total_transcript = round(frac_total_transcript);
    P_count_vec = round(frac_P_count_vec);
    no_mrna = length(type_idx_array); 
    
    if sum(total_transcript) ~= no_mrna
        disp('Indexing Error with total_transcript! Program paused, Press Ctrl-C to exit!');
        pause;
    end
    
    if length(temp) ~= no_mrna
        disp('Indexing Error with temp! Program paused, Press Ctrl-C to exit!');
        pause;
    end
    
    if length(location_array) ~= P_count_vec(1)
        disp('Logical Error: number of ribosomes does not correspond to number of R proteins. Program paused, Press Ctrl-C to exit!');
        pause;
    end
    
    if sum(transition_array<0)>0
        disp('Logical Error: numbers in transition_array are negative. Program paused, Press Ctrl-C to exit!');
        pause;
    end
    
end