function [state_array, location_array, time_out, time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_v5 (state_array, location_array, beta_array_cell, time_in, time_P_cell, P_count_vec, temp, transition_array)
    
    %best version, transition array takes smaller memory and faster, no for
    %loop
    %pool ribosomes have location zero
    %~5 times faster than v1 on test
    
    R = length(state_array);

    if nargin == 7
        no_mrna = length(beta_array_cell);
        transition_array = zeros(no_mrna,R);
        for mrna = 1:no_mrna
            %at first run nothing will be occupied, so can be skipped
            %occupied = unique(state_array(location_array==mrna)); 
            %betas_to_zero = occupied-1; %since if the 2nd place is occupied, moving from first is inhbited
            %beta_array_cell{mrna}(betas_to_zero) = 0;
            transition_array(mrna,:) = beta_array_cell{mrna}(state_array.*(location_array==mrna)+1*(location_array~=mrna)); %otherwise error because e.g. site 20 may not exist on a short mRNA, add one for places that will be zeroed in the next line, otherwise zero indexing error
            transition_array(mrna,:) = transition_array(mrna,:).*(location_array==mrna); %dot multiply so that transition array stays R long for each mRNA
        end
        transition_array = sum(transition_array, 1);
    end
    
    free = sum(location_array==0);
    
    transition_array = [sum(temp)*free, transition_array];
    
    total_out_rate = sum(transition_array);
    
    deltaT = exprnd(1/total_out_rate);
    time_out = time_in + deltaT;
    
    saved_rand = rand;
    
    idx = find(cumsum(transition_array/total_out_rate)>saved_rand, 1, 'first');
    tt = 1:R;
    free_ribos = tt(location_array==0);
    if idx == 1
        saved_rand2 = rand;
        idx2 = find(cumsum(temp/sum(temp))>saved_rand2, 1, 'first');
        ribo_loc = idx2;
        ribo = free_ribos(randperm(length(free_ribos), 1));
    else
        ribo = idx-1;
        ribo_loc = location_array(ribo);
    end
        
    prev_state = state_array(ribo);
    prev_location = location_array(ribo);
    
    %taking care of previous ribos
    if prev_state == 2
        temp(prev_location) = beta_array_cell{prev_location}(1);
    elseif prev_state == 1
        temp(ribo_loc) = 0;
    elseif sum(state_array==(prev_state-1))>0
        if sum(location_array(state_array==(prev_state-1))==prev_location) == 1
            vec = find(state_array==(prev_state-1));
            rib = vec(location_array(state_array==(prev_state-1))==prev_location);
            transition_array(rib+1) = beta_array_cell{prev_location}(prev_state-1);
        end
    end         
    
    %update
    if location_array(ribo)~=0
        if state_array(ribo) == length(beta_array_cell{location_array(ribo)})
            state_array(ribo) = 1;
            location_array(ribo) = 0;
            P_count_vec(ribo_loc) = P_count_vec(ribo_loc)+1;
            time_P_cell{ribo_loc} = [time_P_cell{ribo_loc}, time_out];
        else
            state_array(ribo) = state_array(ribo) + 1;
        end
    else 
        state_array(ribo) = 2;
        location_array(ribo) = ribo_loc;
    end
    
    %taking care of updated transition rate for chosen ribo
    next_state = state_array(ribo);
    next_location = location_array(ribo);
    
    if next_location ~= 0
        if sum(state_array==(next_state+1))>0
            if sum(location_array(state_array==(next_state+1))==next_location) == 0
                transition_array(ribo+1) = beta_array_cell{next_location}(next_state);
            else
                transition_array(ribo+1) = 0;
            end
        else
                transition_array(ribo+1) = beta_array_cell{next_location}(next_state);
        end
    else
        transition_array(ribo+1) = 0;
    end
    transition_array = transition_array(2:end);
end

