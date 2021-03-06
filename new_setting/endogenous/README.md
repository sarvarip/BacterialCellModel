# Explanation of scripts:

#### Endo_WholeCell_final_ns:
The function that is called by all endogenous simulation scripts. It calculates the output at the (t+1)th timestep given the (t)th timestep as input. 
It uses Gillespie-like stochastic simulation. (Gillespie-like because the transition matrix is time-variant: for example, some ribosomes 
will be temporarily restricted from moving when they are blocked by ribosomes ahead - see Figure 2)

#### Endo_WholeCell_final_ns_wrapper:
Simulation script example: It loads a previous endogenous simulation result and continues the endogenous simulation. The n_s (nutrient efficiency) parameter can be tweaked to see how different endogenous simulations evolve. Simply load the relevant script that you want to continue and change n_s so that it corresponds to the n_s used in the simulation that you load. 

#### endo_mRNA_protein_mass_ratio_vs_GR:
Example applications script: produces a rudimentary version of Figure 4B (dependence of mRNA to protein mass ratio on growth rate). 

#### mode_least:
Helper function used by simulation scripts. 

