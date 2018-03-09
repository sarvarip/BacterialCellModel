R0 = 30;
betas = cell(1,4);
betas{1} = ones(1,10); 
betas{1}(1) = 10^-1; 
betas{2} = ones(1,5);
betas{3} = ones(1,10);
betas{4} = ones(1,5);
betas{4}(3) = 0.2;
maxsteps = 40000; 
time = zeros(1,maxsteps+1);
time_P_cell = cell(1,4);
P_count_vec = zeros(1,4);
state_array = ones(1,R0);
location_array = zeros(1,R0);
for timestep=1:maxsteps
    %disp(state_array)
    %disp(location_array)
    [state_array, location_array, time(timestep+1), time_P_cell, P_count_vec] = Gillespie_STS_Prod_Rate_Multi(state_array, location_array, betas, time(timestep), time_P_cell, P_count_vec);
end
time_elapsed = time(end)-time(1001);
transient_P=zeros(1,4);
for i=1:4
    transient_P = sum(time_P_cell{i}<=time(1001)); %we only can about states from the 1001th timstep on (incl), but
    %if a protein was produced at time = 1001, it was the trasition from state
    %1000 to 1001, so we disregard it!
end
    P_ss = P_count_vec-transient_P;
    production_rate = P_ss/time_elapsed;


disp(production_rate);