%% Solution of the steady-state RFM for an mRNA of length 2

clear

lambda = [0.1, 1, 1];
R0 = 10;
soln = solve(@(x) lambda(1)*(R0-x-(lambda(2)*x)/(lambda(3)+x*lambda(2)))*(1-x)-lambda(3)*(lambda(2)*x)/(lambda(3)+x*lambda(2)),'Real', true, 'MaxDegree', 3);
pi_1 = eval(soln(eval(soln)<1&eval(soln)>=0));
pi_2 = (lambda(2)*pi_1)/(lambda(3)+pi_1);
u = R0-pi_1-pi_2;


%% Recursive series solution 

%recursive series solution gives same result

R0 = 10;
syms X
lambda = [0.1, 1, 1]; 
len = length(lambda);

n = sym('n', [len-1,1]); %create symbolic arrays otherwise, cannot convert symbolic to double
d = sym('d', [len-1,1]);
n(1) = lambda(end)*X; %from 3' to 5' now!!
d(1) = lambda(end-1)*(1-X);
for i = 2:(len-1)
    n(i) = lambda(end)*d(i-1)*X;
    d(i) = lambda(end-i)*(d(i-1)-n(i-1));
end
frac = n./d;
F = sum(frac) + X - R0;

soln_series = solve(F, 'Real', true, 'MaxDegree', 3);
pi_2_series = eval(soln_series(eval(soln_series)<1&eval(soln_series)>=0));
X = pi_2_series;
pi_1_series = eval(subs(frac(1)));
u_series = eval(subs(frac(2)));

%% Now try with length 4, 5th order poly


clear

lambda = [0.1, 1, 1, 1]; 
len = length(lambda);
R0 = 10;
syms X

n = sym('n', [len-1,1]); %create symbolic arrays otherwise, cannot convert symbolic to double
d = sym('d', [len-1,1]);
n(1) = lambda(end)*X; %from 3' to 5' now!!
d(1) = lambda(end-1)*(1-X);
for i = 2:(len-1)
    n(i) = lambda(end)*d(i-1)*X;
    d(i) = lambda(end-i)*(d(i-1)-n(i-1));
end
frac = n./d;
F = sum(frac) + X - R0;

soln_series = solve(F, 'MaxDegree', 5);
pi_3_series = vpa(soln_series(vpa(soln_series)<1&vpa(soln_series)>=0));


%Manually solve with toolbox, because Matlab only gives z...sym2poly!
%Use vpa() to numerically approximate the roots!!
% r = roots([1 -9001/750 +14626/375 -9067/200 63001/3000 -10/3]); %last three is possible
% pi_3_series = r(5); %choose last

X = pi_3_series; %0.3820
pi_2_series = eval(subs(frac(1))); %0.618
pi_1_series = eval(subs(frac(2))); %indeed since initiation rate is huge, this is almost 1
u_series = eval(subs(frac(3))); %8

%% Now try with length 5, 6th order poly

clear

lambda = [0.1, 1, 1, 1, 1]; 
len = length(lambda);
R0 = 10;
syms X

n = sym('n', [len-1,1]); %create symbolic arrays otherwise, cannot convert symbolic to double
d = sym('d', [len-1,1]);
n(1) = lambda(end)*X; %from 3' to 5' now!!
d(1) = lambda(end-1)*(1-X);
for i = 2:(len-1)
    n(i) = lambda(end)*d(i-1)*X;
    d(i) = lambda(end-i)*(d(i-1)-n(i-1));
end
frac = n./d;
F = sum(frac) + X - R0;

soln_series = solve(F, 'MaxDegree', 6);
pi_4_series = vpa(soln_series(vpa(soln_series)<1&vpa(soln_series)>=0)); %seems that always the first (smallest soln) is correct

X = pi_4_series; %0.333
pi_3_series = eval(subs(frac(1))); %0.5
pi_2_series = eval(subs(frac(2))); %0.667
pi_1_series = eval(subs(frac(3)));%indeed since initiation rate is huge, this is almost 1
u_series = eval(subs(frac(4))); %7.5

%% length 8

clear

lambda = [0.1, 1, 1, 1, 1, 1, 1, 1]; 
len = length(lambda);
R0 = 10;
syms X

n = sym('n', [len-1,1]); %create symbolic arrays otherwise, cannot convert symbolic to double
d = sym('d', [len-1,1]);
n(1) = lambda(end)*X; %from 3' to 5' now!!
d(1) = lambda(end-1)*(1-X);
for i = 2:(len-1)
    n(i) = lambda(end)*d(i-1)*X;
    d(i) = lambda(end-i)*(d(i-1)-n(i-1));
end
frac = n./d;
F = sum(frac) + X - R0;

soln_series = solve(F, 'MaxDegree', 10);
pi_7_series = vpa(soln_series(vpa(soln_series)<1&vpa(soln_series)>=0)); %seems that always the first (smallest soln) is correct

X = pi_7_series; %0.2831
pi_6_series = eval(subs(frac(1))); %0.3949
pi_5_series = eval(subs(frac(2))); %0.4679
pi_4_series = eval(subs(frac(3))); %0.5321
pi_3_series = eval(subs(frac(4))); %0.6051
pi_2_series = eval(subs(frac(5))); %0.7169
pi_1_series = eval(subs(frac(6)));%indeed since initiation rate is huge, this is almost 1
u_series = eval(subs(frac(7))); 

%% length 20 - 70 solutions! (took ~5 mins)

clear

lambda = [0.1, ones(1,19)]; 
len = length(lambda);
R0 = 25; %changed to 25!
syms X

n = sym('n', [len-1,1]); %create symbolic arrays otherwise, cannot convert symbolic to double
d = sym('d', [len-1,1]);
n(1) = lambda(end)*X; %from 3' to 5' now!!
d(1) = lambda(end-1)*(1-X);
for i = 2:(len-1)
    n(i) = lambda(end)*d(i-1)*X;
    d(i) = lambda(end-i)*(d(i-1)-n(i-1));
end
frac = n./d;
F = sum(frac) + X - R0;

soln_series = solve(F, 'MaxDegree', 21);
pi_19_series = vpa(soln_series(1)); %seems that always the first (smallest soln) is correct

X = pi_19_series; 

pi_n_series = zeros(length(lambda)-1,1);

for j = length(lambda)-1:-1:1
    pi_n_series(j) = eval(subs(frac(length(lambda)-j)));
end

%% Compare this to Gillespie

R0 = 25;
no_sims = 10; %sample stationary distribution instead
maxsteps = 5000; %completely different in state_count(3) in case of 99 than when I set this to a 100 (0.1 vs. 0) -> limit cycle?
ss = 2000;
%coverage = zeros(no_sims,maxsteps+1);
time = zeros(no_sims,maxsteps+1);
state_saved = zeros(maxsteps-ss,R0); %Assuming that stationary distribution is achieved before 100 steps
total_mRNA = 1;
type_mRNA = 1;
betas = cell(1,type_mRNA);
betas{1} = ones(1,20); %3, 4, 5, 8, 20 -> gets better and better as length increases
betas{1}(1) = 0.1; %can increase to infinity, coverage still stays at 0.5, if decrease then happens what expected 
type_idx_array = [1]; %DON'T forget to change this line!
%figure(1)
%hold on
state_count = zeros(length(betas{1}),no_sims); %TAKE TIME INTO ACCOUNT WHEN TAKING AVG!!!
for sim=1:no_sims
    %disp(sim)
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
        if timestep > ss
            state_saved(timestep-ss,:) = state_array;
            %if timestep is 101, then the returned variable is the 102th
            %and the difference (between 102th and 101th) tells you how long we stay in the 101th state
        end
    end
    %plot(time(sim,:),coverage(sim,:))
    for i = 1:length(betas{1})
    state_count(i,sim) = sum(sum(bsxfun(@times,(state_saved==i),(diff(time(sim,(ss+1):end)))'))); %because when first state_array is created time is already 2 and when timestep is 101 (first time it is bigger than 100) time is 102
    end
    state_count(:,sim) = state_count(:,sim)/(R0*(time(sim,end)-time(sim,ss+1)));
end

disp(mean(state_count,2)*R0);

%% Compare protein production rates

R0 = 25;
no_sims = 1; %sample stationary distribution instead
maxsteps = 20000; %completely different in state_count(3) in case of 99 than when I set this to a 100 (0.1 vs. 0) -> limit cycle?
ss = 10000;
%coverage = zeros(no_sims,maxsteps+1);
time = zeros(no_sims,maxsteps+1);
total_mRNA = 1;
type_mRNA = 1;
betas = cell(1,type_mRNA);
betas{1} = ones(1,20); %3, 4, 5, 8, 20 -> gets better and better as length increases
betas{1}(1) = 0.1; %can increase to infinity, coverage still stays at 0.5, if decrease then happens what expected 
type_idx_array = [1]; %DON'T forget to change this line!
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
end
transient_P = sum(time_P_cell{1}<=time(ss+1)); %we only can about states from the 101th timstep on (incl), but 
%if a protein was produced at time = 101, it was the trasition from state
%100 to 101, so we disregard it!
P_ss = P_count_vec-transient_P;
time_elapsed = time(end)-time(ss+1);
production_rate = P_ss/time_elapsed;

disp(production_rate);