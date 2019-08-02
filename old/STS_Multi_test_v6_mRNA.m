%Beta has size of the number of types of mRNA
%Executing lasts a bit longer in this test, but much easier to incorporate
%changing trnacript copies in this version -> useful for whole cell model

R0 = 500; %vs 100; the one that suffers in production is the one with the least initiation rate
total_mRNA = 110;
type_mRNA = 4;
betas = cell(1,type_mRNA);
betas{1} = [0.1, ones(1,9)];
betas{2} = ones(1,5);
betas{3} = ones(1,10);
betas{4} = [ones(1,2), 0.2, ones(1,2)];
type_idx_array = [ones(1,80), 2*ones(1,10), 3, 4*ones(1,19)]; %DON'T forget to change this line!
temp = [0.1*ones(1,80), 1*ones(1,10), 1, 1*ones(1,19)]; %DON'T forget to change this line!
maxsteps = 400000; 
time = zeros(1,maxsteps+1);
time_P_cell = cell(1,total_mRNA);
P_count_vec = zeros(1,total_mRNA);
state_array = ones(1,R0);
location_array = zeros(1,R0);
tic
for timestep=1:maxsteps
    %disp(state_array)
    %disp(location_array)
    if timestep == 1
        [state_array, location_array, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_v6(state_array, location_array, betas, type_idx_array, time(timestep), time_P_cell, P_count_vec, temp);
    else
        [state_array, location_array, time(timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_v6(state_array, location_array, betas, type_idx_array, time(timestep), time_P_cell, P_count_vec, temp, transition_array);
    end
end
toc
time_elapsed = time(end)-time(1001);
transient_P=zeros(1,total_mRNA);
for i=1:total_mRNA
    transient_P(i) = sum(time_P_cell{i}<=time(1001)); %we only can about states from the 1001th timstep on (incl), but
    %if a protein was produced at time = 1001, it was the trasition from state
    %1000 to 1001, so we disregard it!
end
    P_ss = P_count_vec-transient_P;
    production_rate = P_ss/time_elapsed;


disp(production_rate);