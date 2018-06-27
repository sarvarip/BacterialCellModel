figure; 
r = [0.0092, 0.0134, 0.0147, 0.0170, 0.0238, 0.0247, 0.0255];
gr = [0.2035, 0.8815, 1.7030, 2.7323, 4.3379, 4.7687, 4.81];
one = [1 1 1 1 1 1 1];
gr1 = [one; gr]';
regression = gr1\r';
yCalc = regression(2)*gr+regression(1);
scatter(gr,r,200,'filled')
hold on
plot(gr, yCalc, 'LineWidth',4)
set(gca,'fontsize',20)
h = legend('Model results', 'Linear regression');
h.FontSize = 20;
h.Location = 'southeast';
title('Dependence of mRNA/protein ratio on growth rate', 'FontSize', 20);
xlabel('Growth rate (h^-^1)', 'FontSize', 20);
ylabel('mRNA/protein mass ratio', 'FontSize', 20);