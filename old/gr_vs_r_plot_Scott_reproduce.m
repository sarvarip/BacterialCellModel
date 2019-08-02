figure; 
%r = [0.0092, 0.0134, 0.0147, 0.0170, 0.0238, 0.0247, 0.0255];
%r obtained by calculating (total_transcript)*[7500 300 300]' / ((P_count_vec)*[7500 300 300]')
r = [0.009, 0.0129, 0.0147, 0.0174, 0.0221, 0.0248, 0.0254];
std = [2.5757*10^-4, 1.9311*10^-4, 2.7731*10^-4, 8.4069*10^-4, 6.1844*10^-4, 2.9106*10^-4, 5.8602*10^-4];

%to make it with error bars
%res_array = zeros(1,5);
% for i = 2000:2000:10000
% res_array(i/2000) = (transcript_array(i,:))*[7500 300 300]' / (P_count_vec_array(i,:)*[7500 300 300]');
% end
% r = mean(res_array);
% std = std(res_array);

%with more samples - std gets smaller
% vec_rat = transcript_array*[7500 300 300]' ./ (P_count_vec_array*[7500 300 300]');
% r = mean(vec_rat(2:end));
% std = std(vec_rat(2:end));

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
%errorbar(gr,r,std);
set(gca,'fontsize',20)
h = legend('Model results', 'Linear regression');
h.FontSize = 20;
h.Location = 'southeast';
title('Dependence of mRNA/protein ratio on growth rate', 'FontSize', 20);
xlabel('Growth rate (h^-^1)', 'FontSize', 20);
ylabel('mRNA/protein mass ratio', 'FontSize', 20);