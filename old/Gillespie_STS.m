function [state_array, time_out] = Gillespie_STS (state_array, beta_array, time_in)
    occupied = unique(state_array); %will contain 1 (ribosome pool, so take that away)
    if occupied == 1
        occupied = [];
    else
        occupied = occupied(2:end);
    end
    betas_to_zero = occupied-1;
    beta_array(betas_to_zero) = 0;
    total_out_rate = sum(beta_array(state_array));
    deltaT = exprnd(1/total_out_rate);
    time_out = time_in + deltaT;
    
%     comparison_array = 1/total_out_rate * a.*beta_array;

    saved_rand = rand;
    idx = find(cumsum(beta_array(state_array)/total_out_rate)>saved_rand, 1, 'first');
    
%     [val_high,idx_high] = min(comparison_array(comparison_array>saved_rand));
%     diff = val_high-saved_rand;
%     chosen_beta = beta_array(idx_high);
%     increment = chosen_beta/total_out_rate;
%     extra_idx = floor(diff/increment);
%     if extra_idx == 0
%         temp = comparison_array(state_array)==val_high;
%         idx = find(temp,1,'last');
%     else
%         temp = comparison_array(state_array)==val_high;
%         idx = find(temp,extra_idx,'first');
%         idx = idx(end);
%     end

    if state_array(idx) == length(beta_array)
        state_array(idx) = 1;
    else
        state_array(idx) = state_array(idx) + 1;
    end
end

    
    