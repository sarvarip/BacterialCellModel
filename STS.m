clc, clear
close all

%test
% state_array = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %one corresponds to being in ribosome pool
% %in deterministic simulation, -1 coresponds to that
% betas = [1 1 1 1 1 1 1 1 1 10^-5];
% time = 0;
% while time<100
%     [state_array, time] = Gillespie_STS(state_array, betas, time);
%     disp(state_array);
%     disp(time);
%     pause;
% end

%simulation - gonna be very memory intensive -> python, cloud?
betas = ones(1,101); %1001
betas(1) = 10^-4; %can increase to infinity, coverage still stays at 0.5, if decrease then happens what expected 
%betas(2) = 10^-4; %slow codon at begining, correct
%betas(end) = 10^-4; %slow codon at end, correct
slow_loc = 51;
betas(slow_loc) = 100;%0.1 %slow codon massively reduces cell-to-cell variance
no_sims = 3;
maxsteps = 20000;
coverage = zeros(no_sims,maxsteps+1);
coverage_after_slowcod = zeros(no_sims,maxsteps+1);
coverage_before_slowcod = zeros(no_sims,maxsteps+1);
time = zeros(no_sims,maxsteps+1);
figure(1)
hold on
for sim=1:no_sims
    disp(sim)
    state_array = ones(1,10000);
    for timestep=1:maxsteps
        [state_array, time(sim,timestep+1)] = Gillespie_STS(state_array, betas, time(sim,timestep));
        coverage(sim,timestep+1) = sum(state_array>1)/(length(betas)-1);
        coverage_after_slowcod(sim,timestep+1) = sum(state_array>slow_loc)/(length(betas)-slow_loc);
        coverage_before_slowcod(sim,timestep+1) = (sum(state_array>1) - sum(state_array>slow_loc))/(slow_loc-1);
    end
    plot(time(sim,:),coverage(sim,:))
end
%interpolate and average
%need to understand why as initiation rate tends to infinity, coverage is
%still 0.5!!!!!!!
figure(2)
hold on
for sim = 1:no_sims
    plot(time(sim,:),coverage_after_slowcod(sim,:)) %makes sense in terms of deterministic
end

figure(3)
hold on
for sim = 1:no_sims
    plot(time(sim,:),coverage_before_slowcod(sim,:)) %does not really get fully covered until slow codon rate is very slow
end



