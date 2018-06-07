clear

% L = 50;
% M = 3;
% alfaplus = 10^-4;
% alfaminus = 200;
% beta = ones(L+1,1); %everything shifted by one, since there is no idx zero in Matlab
% 
%edge case

tic;

L = 200; %length
M = 1;
alfaplus = 1;
alfaminus = 0;
beta = ones(L+1,1); %first index is free ribosomes
slow_loc = 101;
beta(slow_loc) = 0.1;

tspan = [0 800];
R0 = 1000;
Y_ini = [R0, zeros(1,L+1)];
[t, Y] = ode23(@(t,Y) simulation_Algar(t,Y,L,M,alfaplus,alfaminus,beta), tspan, Y_ini);

toc;

%Plot the density graph - reduce resolution by ~1000 so that matlab can
%actually plot is w/o sinigificant computational issues


%Altogether
time_red = zeros(100,1);
density_red = zeros(100,1);
step = floor(size(Y,1)/100);

figure(1)
hold on

for j = 1:100
    time_red(j) = t(j*step);
    y_red = Y(j*step,:);
    density_red(j) = sum(y_red(2:end))/(L+1);
end
plot(time_red, density_red, 'LineWidth', 3)
%title('Ribosome density on mRNA over time', 'FontSize', 20);
xlabel('Time', 'FontSize', 15);
ylabel('Density', 'FontSize', 15);
%legend({'Average of 10 simulations', 'RFM'}, 'FontSize', 15, 'Location','southeast');
xlim([0 800])

%Before slow codon
time_red = zeros(100,1);
density_red = zeros(100,1);
step = floor(size(Y,1)/100);

figure(2)
hold on

for j = 1:100
    time_red(j) = t(j*step);
    y_red = Y(j*step,:);
    density_red(j) = sum(y_red(2:slow_loc))/((L+1)/2);
end
plot(time_red, density_red, 'LineWidth', 3)
%title('Ribosome density on mRNA over time - Before slow codon', 'FontSize', 20);
xlabel('Time', 'FontSize', 15);
ylabel('Density', 'FontSize', 15);
%legend({'Average of 10 simulations', 'RFM'}, 'FontSize', 15, 'Location','southeast');
xlim([0 800])

%After slow codon
time_red = zeros(100,1);
density_red = zeros(100,1);
step = floor(size(Y,1)/100);

figure(3)
hold on

for j = 1:100
    time_red(j) = t(j*step);
    y_red = Y(j*step,:);
    density_red(j) = sum(y_red(slow_loc+1:end))/((L+1)/2);
end
plot(time_red, density_red, 'LineWidth', 3)


%title('Ribosome density on mRNA over time - After slow codon', 'FontSize', 20);
xlabel('Time', 'FontSize', 15);
ylabel('Density', 'FontSize', 15);
xlim([0 800])
%legend({'Average of 10 simulations', 'RFM'}, 'FontSize', 15, 'Location','southeast');