clear

% L = 50;
% M = 3;
% alfaplus = 10^-4;
% alfaminus = 200;
% beta = ones(L+1,1); %everything shifted by one, since there is no idx zero in Matlab
% 
%edge case

L = 19;
M = 1;
alfaplus = 1;
alfaminus = 0;
beta = ones(L+1,1); %everything shifted by one, since there is no idx zero in Matlab
beta(10) = 0.2;

tspan = [0 300];
Y_ini = [10000, zeros(1,L+1)];
[t, Y] = ode23(@(t,Y) simulation_Algar(t,Y,L,M,alfaplus,alfaminus,beta), tspan, Y_ini);

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
plot(time_red, density_red, 'LineWidth', 2)
title('Ribosome density on mRNA over time', 'FontSize', 20);
xlabel('Time', 'FontSize', 15);
ylabel('Density', 'FontSize', 15);

%Before slow codon
time_red = zeros(100,1);
density_red = zeros(100,1);
step = floor(size(Y,1)/100);

figure(2)
hold on

for j = 1:100
    time_red(j) = t(j*step);
    y_red = Y(j*step,:);
    density_red(j) = sum(y_red(2:11))/((L+1)/2);
end
plot(time_red, density_red, 'LineWidth', 2)
title('Ribosome density on mRNA over time - Before slow codon', 'FontSize', 20);
xlabel('Time', 'FontSize', 15);
ylabel('Density', 'FontSize', 15);

%After slow codon
time_red = zeros(100,1);
density_red = zeros(100,1);
step = floor(size(Y,1)/100);

figure(3)
hold on

for j = 1:100
    time_red(j) = t(j*step);
    y_red = Y(j*step,:);
    density_red(j) = sum(y_red(12:end))/((L+1)/2);
end
plot(time_red, density_red, 'LineWidth', 2)
title('Ribosome density on mRNA over time - After slow codon', 'FontSize', 20);
xlabel('Time', 'FontSize', 15);
ylabel('Density', 'FontSize', 15);
    