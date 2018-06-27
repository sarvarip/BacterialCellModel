codon = 1:19;
occupation_stocell_avg = [0.83, 0.7, 0.64, 0.62, 0.58, 0.56, 0.53, 0.51, 0.50, 0.48, 0.48, 0.47, 0.46, 0.47, 0.44, 0.41, 0.37, 0.33, 0.28];
occupation_rfm = [0.83, 0.69, 0.63, 0.6, 0.57, 0.55, 0.54, 0.53, 0.51, 0.5, 0.49, 0.48, 0.47, 0.46, 0.44, 0.42, 0.39, 0.34, 0.26];

scatter(codon, occupation_stocell_avg, 100, 'filled')
hold on
scatter(codon, occupation_rfm, 'x', 'LineWidth', 3)
set(gca,'fontsize',20)
title('Occupation probabilities of ribosome sites', 'FontSize', 20);
xlabel('Ribosome site from 5 prime end', 'FontSize', 20);
ylabel('Occupation probability', 'FontSize', 20);
h = legend('Model result', 'RFM result');
h.FontSize = 20;
h.Location = 'northeast';