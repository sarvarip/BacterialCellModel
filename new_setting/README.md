# Simulation visualization on YouTube
https://www.youtube.com/watch?v=EPMutCumilk&feature=youtu.be

# Explanation of scripts

### Endogenous results
Please see the folder called endogenous.

### Heterologous results

#### Simulation1 and Simulation2: 
Heterologous simulation scripts that load an endogenous simulation result and then simulate adding a heterologous construct to the bacteria. The two scripts are identical, except for having different index_array and Simulation 2 running twice as many iterations. Simulation 2 includes the cases (specified by index_array) that generally converge slower (see convergence definition), hence the increased number of iterations. We ran these two scripts on HPC servers to get results for the case when the heterologous mRNA has 60 ribosome footprints (as opposed to 30, see below). Note that depending on the available nodes, we changed index_array, but left the rest of the script unchanged. We used these scripts to get a converged output for each of the cases (9 endogenous, 9 heterologous with slow codon at 5th position and 9 heterologous with slow codon at 55th position). 

#### Whole_cell_final_par:
The function that is called by all simulation scripts. It calculates the output at the (t+1)th timestep given the (t)th timestep as input. It uses Gillespie-like stochastic simulation. (Gillespie-like because the transition matrix is time-variant: for example, some ribosomes will be temporarily restricted from moving when they are blocked by ribosomes ahead - see Figure 2)

#### Whole_cell_final_par_wrapper, Whole_cell_final_par_wrapper_003, Whole_cell_final_par_wrapper_0005:
These are heterologous simulation scripts, very similar to Simulation1 and Simulation 2 above. The difference is that these scripts consider the case when the heterologous mRNA has 30 ribosome places (as opposed to 60, see above). In addition, Whole_cell_final_par_newsetting considers slow codons with relative speeds of 0.03 and 0.005. Whole_cell_final_par_newsetting has been broken down to two scripts that were run on two different servers. Whole_cell_final_par_newsetting_003 considers only the 0.03 slow codon speed case, whereas Whole_cell_final_par_newsetting_0005 considers only the 0.005 slow codon speed case. The results are clearer in the latter case as the slower codon speed emhasizes the ribosome sequestration and burden. 

#### mode_least and parsave_final:
Helper functions that are called by the simulation scripts. 

# Explanation of variables

#### State_array: 
Has length equal to the number of ribosomes; describes which ribosome footprint the ribosome is translating at the moment; if ribosome is free, the respective number is one, if ribosome is at the beginning of the transcript, its state is two and if it is at the nth position, its state is n+1. Note: this does not cause an issue, because if transcript has length n, beta array (see below) of that transcript is n+1 long, so protein is produced when the state is equal to the length of the corresponding beta array!

#### Location_array: 
Has length equal to the number of ribosomes; describes which mRNA the ribosome is on, if the ribosome is free, the respective number is zero

#### Time_P_cell: 
Matlab cellular array, which has length equal to the number of types of mRNAs; the respective arrays store the times when a protein type was made from a given mRNA type

#### P_count_vec: 
Array, has length equal to the number of protein types; the respective elements store the total number of proteins of the corresponding type

#### Betas: 
Matlab cellular array, which has length equal to the number of types of mRNAs; the respective arrays store the local elongation speed for each ribosome footprint on a given mRNA starting with initiation speed and hence have the length of the length of the mRNA plus one (because of initiation speed) Note: initiation rate should be reflected both in beta_array_cell and in temp!

#### Transition_array – after its creation, which happens when number of input arguments to Whole_cell_final_par is 12 (as opposed to 13) because transition array is not passed in, so we need to define it (line 19): 
Matrix (2D array), which has number of rows equal to the number of mRNAs and number of columns equal to the number of ribosomes; each row stores the current elongation speeds associated with the ribosomes that are on the respective mRNA, otherwise the element is zero

#### Transition_array – after merge and sum: 
1D array (as opposed to it being initially a matrix - see above), which has length of the number of ribosomes plus eight. The first element (created using an array called temp) corresponds to the total transition rate of the free ribosomes to any of the mRNAs, whereas the following elements (each element corresponding to one ribosome) correspond to the transition rate of a ribosome to the next position on the mRNA; if the ribosome is free, the element is zero. The seven last elements correspond to R transcription rate, E transcription rate, Q transcription rate, H transcription rate, degradation rate, nutrient transport rate and nutrient metabolism rate, respectively. 

#### Temp: 
Array, which has length equal to the number of mRNAs. Each element corresponds to the current initiation rate of the mRNAs taking into account that if the RBS site is occupied the current initiation rate is zero

#### Type_idx_array: 
Array, which has length equal to the number of mRNAs. Each element denotes which type the corresponding mRNA belongs to. 

# Simple applications

#### Get growth rate
`avg_inst_growth_rate`

#### Get heterologous protein production
`production_rate(4)`

#### Get ribosome density:
`state_a = state_array(state_array>1);
location_a = location_array(location_array>0);
hist(state_a(type_idx_array(location_a)==4)-1,30)`

#### Heterologous mass fraction calculation (mean over some time range after convergence is achieved):
`100 * (mean(P_count_vec_array(36000:end,4))*300) / ( mean(P_count_vec_array(36000:end,1))*7500 + mean(P_count_vec_array(36000:end,2))*300 + mean(P_count_vec_array(36000:end,3))*300 + mean(P_count_vec_array(36000:end,4))*300 )`

#### Heterologous mass fraction standard deviation:
`100 * (mean(P_count_vec_array(36000:end,4))*300) / ( mean(P_count_vec_array(36000:end,1))*7500 + mean(P_count_vec_array(36000:end,2))*300 + mean(P_count_vec_array(36000:end,3))*300 + mean(P_count_vec_array(36000:end,4))*300 )`

#### Check for convergence:
`th = mean(P_count_vec_array(0.9*end:end,4))
upper = 1.01*th
lower = 0.99 * th
plot(time_ss, P_count_vec_array(:,4))
hold on
plot(time_ss, upper*ones(length(time_ss),1))
plot(time_ss, lower*ones(length(time_ss),1))
plot(0.9*time_ss(end)*ones(length(time_ss),1), linspace(1,P_count_vec_array(end,4),length(time_ss)))`






