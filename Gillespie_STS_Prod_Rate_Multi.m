function [state_array, location_array, time_out, time_P_new_cell, P_count_new_cell] = Gillespie_STS_Prod_Rate_Multi (state_array, location_array, beta_array_cell, time_in, time_P_cell, P_count_cell)
    
    %pool ribosomes have location zero

    no_mrna = length(beta_array_cell);
    mrna_length_vec = zeros(1,no_mrna);

    for mrna = 1:no_mrna
        occupied = unique(state_array(location_array==mrna)); 
        occupied = occupied(2:end);
        betas_to_zero = occupied-1; %since if the 2nd place is occupied, moving from first is inhbited
        beta_array_cell{mrna}(betas_to_zero) = 0;
        beta_array_cell{mrna} = beta_array_cell{mrna}(state_array(location_array==mrna | location_array==0));
        mrna_length_vec(mrna) = length(beta_array_cell{mrna});
    end
    
    transition_array = cell2mat(beta_array_cell);
    total_out_rate = sum(transition_array);
    
    deltaT = exprnd(1/total_out_rate);
    time_out = time_in + deltaT;
    
    saved_rand = rand;
    
    idx = find(cumsum(transition_array/total_out_rate)>saved_rand, 1, 'first');
    
    if state_array(idx) == length(beta_array_cell{location_array(idx)})
        state_array(idx) = 1;
        location_array(idx) = 0;
        ribo_loc = find(cumsum(mrna_length_vec)==idx, 1, 'first'); %equal since ribo is on last site of mRNA
        P_count_new_cell{ribo_loc} = P_count_cell{ribo_loc}+1;
        time_P_new_cell{ribo_loc} = [time_P_cell{ribo_loc}, time_out];
    elseif location_array(idx)==0
        ribo_loc = find(cumsum(mrna_length_vec)>idx, 1, 'first'); %definitely bigger, if it was equal, then ribosome would on last site on mRNA
        state_array(idx) = 2;
        location_array(idx) = ribo_loc;
        time_P_new_cell = time_P_cell;
        P_count_new_cell = P_count_cell;
    else
        state_array(idx) = state_array(idx) + 1;
        time_P_new_cell = time_P_cell;
        P_count_new_cell = P_count_cell;
    end
end
