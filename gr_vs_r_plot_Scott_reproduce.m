figure; 
r = [0.0092, 0.0134, 0.0147, 0.0170, 0.0238, 0.0247, 0.0255];
%r obtained by calculating (total_transcript)*[7500 300 300]' / ((P_count_vec)[7500 300 300]')
scaling_factor = 330*3/118.9; %Average weight of condon triplet in E.Coli divided by 
%the average weight of an amino acid in E.Coli
r = r*scaling_factor;
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