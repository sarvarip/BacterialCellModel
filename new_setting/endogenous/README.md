# Explanation of scripts:

#### Gillespie_STS_Prod_Rate_Multi_WholeCell_2_final_ns_10:
The function that is called by all endogenous simulation scripts. It calculates the output at the (t+1)th timestep given the (t)th timestep as input. 
It uses Gillespie-like stochastic simulation. (Gillespie-like because the transition matrix is time-variant: for example, some ribosomes 
will be temporarily restricted from moving when they are blocked by ribosomes ahead - see Figure 2)

#### WholeCellTest_gr_newss_ns_10_transcript_array:
Simulation script example: It loads a previous endogenous simulation result and continues the endogenous simulation. Even though the name has ns_10 in it, the n_s (nutrient efficiency) parameter can be tweaked to see how different endogenous simulations evolve. Simply load the relevant script that you want to continue and change n_s so that it corresponds to the n_s used in the simulation that you load. 

#### endo_mRNA_protein_mass_ratio_vs_GR:
Example applications script: produces a rudimentary version of Figure 4B (dependence of mRNA to protein mass ratio on growth rate). 

#### mode_least:
Helper function used by simulation scripts. 

