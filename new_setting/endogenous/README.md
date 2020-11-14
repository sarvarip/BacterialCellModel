# Explanation of scripts:

#### Gillespie_STS_Prod_Rate_Multi_WholeCell_2_final_ns_10:
The function that is called by all simulation scripts. It calculates the output at the (t+1)th timestep given the (t)th timestep as input. 
It uses Gillespie-like stochastic simulation. (Gillespie-like because the transition matrix is time-variant: for example, some ribosomes will be temporarily restricted from moving when they are blocked by ribosomes ahead - see Figure 2)
