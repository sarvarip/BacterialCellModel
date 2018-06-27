figure; 
r = [0.0092, 0.0134, 0.0147, 0.0170, 0.0238, 0.0247, 0.0255];
gr = [0.2035, 0.8815, 1.7030, 2.7323, 4.3379, 4.7687, 4.81];
one = [1 1 1 1 1 1 1];
gr1 = [one; gr]';
regression = gr1\r';
yCalc = regression(2)*gr+regression(1);
scatter(gr,r,'LineWidth',2)
hold on
plot(gr, yCalc, 'LineWidth',2)
h = legend('Model results', 'Linear regression');
h.FontSize = 15;
h.Location = 'southeast';
title('Dependence of mRNA/protein ratio on growth rate', 'FontSize', 20);
xlabel('Growth rate', 'FontSize', 15);
ylabel('mRNA/protein ratio', 'FontSize', 15);
