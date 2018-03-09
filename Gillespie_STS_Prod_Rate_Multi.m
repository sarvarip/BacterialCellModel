function [state_array, location_array, time_out, time_P_cell, P_count_vec] = Gillespie_STS_Prod_Rate_Multi (state_array, location_array, beta_array_cell, time_in, time_P_cell, P_count_vec)
    
    %pool ribosomes have location zero

    no_mrna = length(beta_array_cell);
    transition_array_cell = cell(1,no_mrna);
    R = length(state_array);
    %mrna_length_vec = zeros(1,no_mrna);

    for mrna = 1:no_mrna
        occupied = unique(state_array(location_array==mrna)); 
        betas_to_zero = occupied-1; %since if the 2nd place is occupied, moving from first is inhbited
        beta_array_cell{mrna}(betas_to_zero) = 0;
        transition_array_cell{mrna} = beta_array_cell{mrna}(state_array.*(location_array==mrna | location_array==0)+1*(location_array~=mrna & location_array~=0)); %otherwise error because e.g. site 20 may not exist on a short mRNA, add one for places that will be zeroed in the next line, otherwise zero indexing error
        transition_array_cell{mrna} = transition_array_cell{mrna}.*(location_array==mrna | location_array==0); %dot multiply so that transition array stays R long for each mRNA
        %mrna_length_vec(mrna) = length(transition_array_cell{mrna});
    end
    
    transition_array = cell2mat(transition_array_cell);
    total_out_rate = sum(transition_array);
    
    deltaT = exprnd(1/total_out_rate);
    time_out = time_in + deltaT;
    
    saved_rand = rand;
    
    idx = find(cumsum(transition_array/total_out_rate)>saved_rand, 1, 'first');
    %ribo_loc = find(cumsum(mrna_length_vec)>=idx, 1, 'first');
    %site = idx-sum(mrna_length_vec(1:ribo_loc-1));
    ribo = rem(idx,R);
    if ribo == 0
        ribo = R;
    end
    if idx/R == floor(idx/R)
        ribo_loc = idx/R;
    else
        ribo_loc = floor(idx/R)+1;
    end
    
    
    if location_array(ribo)~=0
        if state_array(ribo) == length(beta_array_cell{location_array(ribo)})
            state_array(ribo) = 1;
            location_array(ribo) = 0;
            %disp(ribo_loc)
            P_count_vec(ribo_loc) = P_count_vec(ribo_loc)+1;
            time_P_cell{ribo_loc} = [time_P_cell{ribo_loc}, time_out];
        else
            state_array(ribo) = state_array(ribo) + 1;
        end
    else 
        state_array(ribo) = 2;
        location_array(ribo) = ribo_loc;
    end
end
