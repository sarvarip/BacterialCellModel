avg_ratio = zeros(7,1);
std_ratio = zeros(7,1);
load('FYP_22_12_endo_gr_newss_ns_10_transcript_arr')
avg_ratio(1) = mean(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
std_ratio(1) = std(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
load('FYP_22_12_endo_gr_newss_ns_50_transcript_arr')
avg_ratio(2) = mean(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
std_ratio(2) = std(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
load('FYP_22_12_endo_gr_newss_ns_100_transcript_arr')
avg_ratio(3) = mean(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
std_ratio(3) = std(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
load('FYP_22_12_endo_gr_newss_ns_200_transcript_arr')
avg_ratio(4) = mean(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
std_ratio(4) = std(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
load('FYP_22_12_endo_gr_newss_ns_400_transcript_arr')
avg_ratio(5) = mean(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
std_ratio(5) = std(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
load('FYP_22_12_endo_gr_newss_ns_500_transcript_arr')
avg_ratio(6) = mean(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
std_ratio(6) = std(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
load('FYP_22_12_endo_gr_newss_ns_600_transcript_arr')
avg_ratio(7) = mean(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;
std_ratio(7) = std(sum(transcript_array(2:end,:)*[7500 300 300]')/sum(P_count_vec_array(2:end,:)*[7500 300 300]'))*990/118.9;

figure; 
r = avg_ratio;
gr = [0.2035, 0.8815, 1.7030, 2.7323, 4.3379, 4.7687, 4.81];
one = [1 1 1 1 1 1 1];
gr1 = [one; gr]';
regression = gr1\r;
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
