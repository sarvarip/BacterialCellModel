function parsave_final_transition_arr(str, time_ss, state_array, location_array, type_idx_array, temp, energy, energy_array, S_i, total_transcript, total_inst_gr_array, P_count_vec_array, total_P, transient_P, P_ss, time_elapsed, production_rate, growth_rate, avg_inst_growth_rate, std_inst_growth_rate, transition_array)
save(str, 'time_ss', 'state_array', 'location_array', 'type_idx_array', 'temp', 'energy', 'energy_array', 'S_i', 'total_transcript', 'total_inst_gr_array', 'P_count_vec_array', 'total_P', 'transient_P', 'P_ss', 'time_elapsed', 'production_rate', 'growth_rate', 'avg_inst_growth_rate', 'std_inst_growth_rate', 'transition_array')
end

%same as parsave_final_run2