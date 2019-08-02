%Scott paper growth rate vs. mass fraction of unnecessary protein
%order 0.005 case: J,K,L,M,N,P,S,T,U,V,W,Y
%then 0.03 case: A,B,C,D,E,G,J,K,L,M,P,S,T,U,V,W,Y

gr_array = [1.661, 1.514, 0.761, 1.242, 0.557, 0.187, 1.606, 1.426, 0.738, 1.065, 0.569, 0.126, ...
    1.731, 1.599, 1.003, 1.471, 1.048, 1.002, 1.648, 1.530, 1.021, 1.499, 0.786, 1.667, 1.59, ...
    0.955, 1.467, 0.977, 0.665];

percentage = [0.626, 6.000, 42.29, 4.972, 41.38, 35.51, 0.627, 6.303, 42.46, 6.376, 41.04, 41.10, ...
    0.635, 5.942, 38.02, 5.765, 37.48, 43.21, 0.615, 5.827, 38.82, 5.576, 37.40, 5.702, ...
    5.983, 38.29, 5.702, 38.88, 40.46];

%1: no slow, 2: 0.03 3: 0.005
codon = [0.005*ones(1,12), ones(1,6), 0.03*ones(1,11)];

one = ones(1,29);
percentage1 = [one; percentage]';
regression = percentage1\gr_array';
yCalc = regression(2)*percentage+regression(1);
s = gscatter(percentage, gr_array, codon);
s(1).MarkerSize = 25;
s(2).MarkerSize = 25;
s(3).MarkerSize = 25;
hold on
plot(percentage, yCalc, 'LineWidth',4)
set(gca,'fontsize',20)
h = legend('Very slow codon in transcript', 'Slow codon in transcript', 'No slow codon in transcript', 'Linear regression of all data');
h.FontSize = 20;
h.Location = 'southwest';
title('Dependence of growth rate on heterologous protein expression', 'FontSize', 20);
xlabel('Mass fraction of heterologous protein', 'FontSize', 20);
ylabel('Growth rate (h^-^1)', 'FontSize', 20);
    