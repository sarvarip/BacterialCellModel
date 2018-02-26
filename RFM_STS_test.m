%% Solution of the steady-state RFM for an mRNA of length 2

clear

% lambda_arr = [10^3, 10^-3, 1, 1];
% lambda_1_arr = [1, 1, 10^3, 1];
% lambda_2_arr = [1, 1, 1, 10^3];
% u = zeros(1,4);
% pi_1 = zeros(1,4);
% pi_2 = zeros(1,4);
% R0 = 100;
% 
% for cases = 1:4
%     lambda = lambda_arr(cases);
%     lambda_1 = lambda_1_arr(cases);
%     lambda_2 = lambda_2_arr(cases);
%     soln = solve(@(x) lambda*(R0-x-(lambda_1*x)/(lambda_2+x*lambda_1))*(1-x)-lambda_2*(lambda_1*x)/(lambda_2+x*lambda_1),'Real', true, 'MaxDegree', 3);
%     pi_1(cases) = soln(eval(soln)<1&eval(soln)>=0);
%     pi_2(cases) = (lambda_1*pi_1(cases))/(lambda_2+pi_1(cases));
%     u(cases) = 1-pi_1(cases)-pi_2(cases);
% end

lambda = [10^3, 1, 1];
R0 = 10;
soln = solve(@(x) lambda(1)*(R0-x-(lambda(2)*x)/(lambda(3)+x*lambda(2)))*(1-x)-lambda(3)*(lambda(2)*x)/(lambda(3)+x*lambda(2)),'Real', true, 'MaxDegree', 3);
pi_1 = eval(soln(eval(soln)<1&eval(soln)>=0));
pi_2 = (lambda(2)*pi_1)/(lambda(3)+pi_1);
u = R0-pi_1-pi_2;


%% Recursive series solution 

%recursive series solution gives same result

R0 = 10;
syms X
lambda = [10^3, 1, 1]; 
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

lambda = [10^3, 1, 1, 1]; 
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

%% Gillespie for comparison

%put graphs on same plot!

%simulation - gonna be very memory intensive -> python, cloud?
R0 = 10;
betas = ones(1,8); %3, 4, 5, 8 -> gets better and better as length increases
betas(1) = 10^3; %can increase to infinity, coverage still stays at 0.5, if decrease then happens what expected 
no_sims = 1; %sample stationary distribution instead
maxsteps = 10000; %completely different in state_count(3) in case of 99 than when I set this to a 100 (0.1 vs. 0) -> limit cycle?
%coverage = zeros(no_sims,maxsteps+1);
time = zeros(no_sims,maxsteps+1);
state_saved = zeros(maxsteps-100,R0); %Assuming that stationary distribution is achieved before 100 steps
%figure(1)
%hold on
for sim=1:no_sims
    %disp(sim)
    state_array = ones(1,R0);
    for timestep=1:maxsteps(rem(sim,10))
        [state_array, time(sim,timestep+1)] = Gillespie_STS(state_array, betas, time(sim,timestep));
        %coverage(sim,timestep+1) = sum(state_array>1)/(length(betas)-1);
        %disp(state_array)
        if timestep > 100
            state_saved(timestep-100,:) = state_array;
        end
    end
    %plot(time(sim,:),coverage(sim,:))
end

state_count = zeros(length(betas),1); %TAKE TIME INTO ACCOUNT WHEN TAKING AVG!!!
for i = 1:length(betas)
    state_count(i) = sum(sum(bsxfun(@times,(state_saved==i),(diff(time(101:end)))'))); %because when first state_array is created time is already 2 and when timestep is 101 (first time it is bigger than 100) time is 102
end

state_count = (maxsteps-100)*state_count/((maxsteps-100)*R0*(time(end)-time(101)));

%title('Modified Gillespie simulation: Ribosome density', 'FontSize', 20)
%xlabel('Time', 'FontSize', 15)
%ylabel('Ribosome density', 'FontSize', 15)

%% Now try with length 5, 6th order poly

clear

lambda = [10^3, 1, 1, 1, 1]; 
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

lambda = [10^3, 1, 1, 1, 1, 1, 1, 1]; 
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

lambda = [10^3, ones(1,19)]; 
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
betas = ones(1,20); %3, 4, 5, 8, 20 -> gets better and better as length increases
betas(1) = 10^3; %can increase to infinity, coverage still stays at 0.5, if decrease then happens what expected 
no_sims = 1; %sample stationary distribution instead
maxsteps = 10000; %completely different in state_count(3) in case of 99 than when I set this to a 100 (0.1 vs. 0) -> limit cycle?
%coverage = zeros(no_sims,maxsteps+1);
time = zeros(no_sims,maxsteps+1);
state_saved = zeros(maxsteps-100,R0); %Assuming that stationary distribution is achieved before 100 steps
%figure(1)
%hold on
for sim=1:no_sims
    %disp(sim)
    state_array = ones(1,R0);
    for timestep=1:maxsteps(rem(sim,10))
        [state_array, time(sim,timestep+1)] = Gillespie_STS(state_array, betas, time(sim,timestep));
        %coverage(sim,timestep+1) = sum(state_array>1)/(length(betas)-1);
        %disp(state_array)
        if timestep > 100
            state_saved(timestep-100,:) = state_array;
            %if timestep is 101, then the returned variable is the 102th
            %and the difference (between 102th and 101th) tells you how long we stay in the 101th state
        end
    end
    %plot(time(sim,:),coverage(sim,:))
end

state_count = zeros(length(betas),1); %TAKE TIME INTO ACCOUNT WHEN TAKING AVG!!!
for i = 1:length(betas)
    state_count(i) = sum(sum(bsxfun(@times,(state_saved==i),(diff(time(101:end)))'))); %because when first state_array is created time is already 2 and when timestep is 101 (first time it is bigger than 100) time is 102
end

state_count = state_count/(R0*(time(end)-time(101)));

disp(state_count*R0);

%% Compare protein production rates

R0 = 25;
betas = ones(1,20); %3, 4, 5, 8, 20 -> gets better and better as length increases
betas(1) = 10^3; %can increase to infinity, coverage still stays at 0.5, if decrease then happens what expected 
maxsteps = 10000; %completely different in state_count(3) in case of 99 than when I set this to a 100 (0.1 vs. 0) -> limit cycle?
time = zeros(1,maxsteps+1);
time_P = [];
P_count = 0;
state_array = ones(1,R0);
for timestep=1:maxsteps
    [state_array, time(timestep+1), time_P, P_count] = Gillespie_STS_Prod_Rate(state_array, betas, time(timestep), time_P, P_count);
end
transient_P = sum(time_P<=time(101)); %we only can about states from the 101th timstep on (incl), but 
%if a protein was produced at time = 101, it was the trasition from state
%100 to 101, so we disregard it!
P_ss = P_count-transient_P;
time_elapsed = time(end)-time(101);
production_rate = P_ss/time_elapsed;

disp(production_rate);