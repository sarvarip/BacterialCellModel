%Beta has size of the number of types of mRNA
%Executing lasts a bit longer in this test, but much easier to incorporate
%changing transcript copies in this version -> useful for whole cell model

close all
clc, clear

%% Multiple trajectories

tic

R0 = 1000; %vs 100; the one that suffers in production is the one with the least initiation rate
total_mRNA = 1;
type_mRNA = 1;
betas = cell(1,type_mRNA);
betas{1} = [ones(1,99), 0.1, ones(1,100)];
slow_loc = 100;
type_idx_array = [1]; %DON'T forget to change this line!
maxsteps = 15000; 
no_sims = 10;
coverage = zeros(no_sims,maxsteps+1);
coverage_after_slowcod = zeros(no_sims,maxsteps+1);
coverage_before_slowcod = zeros(no_sims,maxsteps+1);
time = zeros(no_sims,maxsteps+1);

for sim=1:no_sims
    disp(sim)
    time_P_cell = cell(1,total_mRNA);
    P_count_vec = zeros(1,total_mRNA);
    state_array = ones(1,R0);
    location_array = zeros(1,R0);
    temp = [1]; %DON'T forget to change this line!
    for timestep=1:maxsteps
        if timestep == 1
            [state_array, location_array, time(sim,timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_v6(state_array, location_array, betas, type_idx_array, time(sim,timestep), time_P_cell, P_count_vec, temp);
        else
            [state_array, location_array, time(sim,timestep+1), time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_v6(state_array, location_array, betas, type_idx_array, time(sim,timestep), time_P_cell, P_count_vec, temp, transition_array);
        end
        coverage(sim,timestep+1) = sum(state_array>1)/(length(betas{1})-1);
        coverage_after_slowcod(sim,timestep+1) = sum(state_array>slow_loc)/(length(betas{1})-slow_loc);
        coverage_before_slowcod(sim,timestep+1) = (sum(state_array>1) - sum(state_array>slow_loc))/(slow_loc-1);
    end
end

toc

figure(1)
hold on
for sim = 1:no_sims
    plot(time(sim,:),coverage(sim,:))
end

%Avg
%plot(time(sim,:), mean(coverage, 1), 'b', 'LineWidth', 3); %You cannot 
%really take the average of stochastic process plots..each element in
%matrix corresponds to different times

figure(2)
hold on
for sim = 1:no_sims
    plot(time(sim,:),coverage_before_slowcod(sim,:)) %does not really get fully covered until slow codon rate is very slow
end

figure(3)
hold on
for sim = 1:no_sims
    plot(time(sim,:),coverage_after_slowcod(sim,:)) %makes sense in terms of deterministic
end

%% Avg plots

% tic
% 
% R0 = 1000; %vs 100; the one that suffers in production is the one with the least initiation rate
% total_mRNA = 1;
% type_mRNA = 1;
% betas = cell(1,type_mRNA);
% betas{1} = [ones(1,99), 0.1, ones(1,100)];
% slow_loc = 100;
% type_idx_array = [1]; %DON'T forget to change this line!
% no_sims = 10;
% 
% time_res=0:1:1000;
% coverage_cell = cell(length(time_res),1);
% coverage_after_slowcod_cell = cell(length(time_res),1);
% coverage_before_slowcod_cell = cell(length(time_res),1);
% 
% for sim=1:no_sims
%     time = 0;
%     time_P_cell = cell(1,total_mRNA);
%     P_count_vec = zeros(1,total_mRNA);
%     state_array = ones(1,R0);
%     location_array = zeros(1,R0);
%     temp = [1]; %DON'T forget to change this line!
%     [state_array, location_array, time, time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_v6(state_array, location_array, betas, type_idx_array, time, time_P_cell, P_count_vec, temp);
%     idx = find(time_res < time, 1, 'last');
%     coverage_cell{idx}(end+1) = sum(state_array>1)/(length(betas{1})-1);
%     coverage_after_slowcod_cell{idx}(end+1) = sum(state_array>slow_loc)/(length(betas{1})-slow_loc);
%     coverage_before_slowcod_cell{idx}(end+1) = (sum(state_array>1) - sum(state_array>slow_loc))/(slow_loc-1);
%     while time < time_res(end)
%         [state_array, location_array, time, time_P_cell, P_count_vec, temp, transition_array] = Gillespie_STS_Prod_Rate_Multi_v6(state_array, location_array, betas, type_idx_array, time, time_P_cell, P_count_vec, temp, transition_array);
%         idx = find(time_res < time, 1, 'last'); 
%         coverage_cell{idx}(end+1) = sum(state_array>1)/(length(betas{1})-1);
%         coverage_after_slowcod_cell{idx}(end+1) = sum(state_array>slow_loc)/(length(betas{1})-slow_loc);
%         coverage_before_slowcod_cell{idx}(end+1) = (sum(state_array>1) - sum(state_array>slow_loc))/(slow_loc-1);
%     end
% end
% 
% toc
% 
% figure(1)
% hold on
% coverage_avg = zeros(1,length(coverage_cell));
% for t = 1:length(coverage_cell)
%     coverage_avg(t) = mean(coverage_cell{t});
% end
% scatter(time_res,coverage_avg, 'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1)
% 
% figure(2)
% hold on
% coverage_before_slowcod_avg = zeros(1,length(coverage_before_slowcod_cell));
% for t = 1:length(coverage_before_slowcod_cell)
%     coverage_before_slowcod_avg(t) = mean(coverage_before_slowcod_cell{t});
% end
% scatter(time_res,coverage_before_slowcod_avg, 'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1)
% 
% figure(3)
% hold on
% coverage_after_slowcod_avg = zeros(1,length(coverage_after_slowcod_cell));
% for t = 1:length(coverage_after_slowcod_cell)
%     coverage_after_slowcod_avg(t) = mean(coverage_after_slowcod_cell{t});
% end
% scatter(time_res,coverage_after_slowcod_avg, 'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1)


% time_elapsed = time(end)-time(1001);
% transient_P=zeros(1,total_mRNA);
% for i=1:total_mRNA
%     transient_P(i) = sum(time_P_cell{i}<=time(1001)); %we only can about states from the 1001th timstep on (incl), but
%     %if a protein was produced at time = 1001, it was the trasition from state
%     %1000 to 1001, so we disregard it!
% end
%     P_ss = P_count_vec-transient_P;
%     production_rate = P_ss/time_elapsed;
% 
% 
% disp(production_rate);