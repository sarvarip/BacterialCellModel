function [state_array, time_out, time_P_new, P_count_new] = Gillespie_STS_Prod_Rate (state_array, beta_array, time_in, time_P, P_count)
    occupied = unique(state_array); %will contain 1 (ribosome pool, so take that away)
    if occupied == 1
        occupied = [];
    else
        occupied = occupied(2:end);
    end
    betas_to_zero = occupied-1; %since if the 2nd place is occupied, moving from first is inhbited
    beta_array(betas_to_zero) = 0;
    total_out_rate = sum(beta_array(state_array));
    deltaT = exprnd(1/total_out_rate);
    time_out = time_in + deltaT;
    
    saved_rand = rand;
    idx = find(cumsum(beta_array(state_array)/total_out_rate)>saved_rand, 1, 'first');
    
    if state_array(idx) == length(beta_array)
        state_array(idx) = 1;
        P_count_new = P_count+1;
        time_P_new = [time_P, time_out];
    else
        state_array(idx) = state_array(idx) + 1;
        time_P_new = time_P;
        P_count_new = P_count;
    end
end

    
    