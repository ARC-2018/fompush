% load('SimulationResults/GeneralTrajectoryLearning/04_16_2017__08_44_21/cost_data_for_simnumber_1000_straight_line.mat')
% load('SimulationResults/GeneralTrajectoryLearning/04_16_2017__08_44_21/std_cost_data_for_simnumber_1000_straight_line.mat')
% 
% load('SimulationResults/GeneralTrajectoryLearning/04_16_2017__05_13_09_cost_data_for_simnumber_50_straight_line.mat')
% load('SimulationResults/GeneralTrajectoryLearning/04_16_2017__05_13_09_std_cost_data_for_simnumber_50_straight_line.mat')

load('SimulationResults/LinePusherLearning/04_21_2017__04_08_44/cost_data_for_simnumber_1000_straight_line.mat')
load('SimulationResults/LinePusherLearning/04_21_2017__04_08_44/std_cost_data_for_simnumber_1000_straight_line.mat')

a = costs_matrix(1,:);
c = costs_matrix(3,:);
d = costs_matrix(4,:);

wid = 0.1;
cou = 20;

figure
h1 = histogram(a, cou);
hold on
h1.BinWidth = wid;
h3 = histogram(c, cou);
h3.BinWidth = wid;
h4 = histogram(d, cou);
h4.BinWidth = wid;
title('Histogram Over Modes Prior Standardization', 'Interpreter', 'LaTex','FontSize', 20)
xlabel(strcat('Value of $C_i(\bar{x_j})$ over all possible i. Box width: ', num2str(wid), '. Box count: ', num2str(cou)), 'Interpreter', 'LaTex','FontSize', 20)
ylabel(strcat('Samples in each value range from ', num2str(size(costs_matrix, 1)),' experiments.'), 'Interpreter', 'LaTex','FontSize', 20)
legend('j = 1', 'j = 3', 'j = 4')

a = standardized_cost_matrix(1,:);
c = standardized_cost_matrix(3,:);
d = standardized_cost_matrix(4,:);

figure
h1 = histogram(a, cou);
hold on
h1.BinWidth = wid;
h3 = histogram(c, cou);
h3.BinWidth = wid;
h4 = histogram(d, cou);
h4.BinWidth = wid;
title('Histogram Over Modes After Standardization', 'Interpreter', 'LaTex','FontSize', 20)
xlabel(strcat('Value of $\bar{C_i}(\bar{x_j})$ over all possible i. Box width: ', num2str(wid), '. Box count: ', num2str(cou)), 'Interpreter', 'LaTex','FontSize', 20)
ylabel(strcat('Samples in each value range from ', num2str(size(costs_matrix, 1)),' experiments.'), 'Interpreter', 'LaTex','FontSize', 20)
legend('j = 1', 'j = 3', 'j = 4')

a = standardized_cost_matrix(:,82);
c = standardized_cost_matrix(:,163);
figure
hold on
h1 = plot(a);
h3 = plot(c);
title('Cost Per Experiments After Standardization', 'Interpreter', 'LaTex','FontSize', 20)
xlabel('Experiment i', 'Interpreter', 'LaTex','FontSize', 20)
ylabel('Value of $\bar{C_i}(\bar{x_j})$', 'Interpreter', 'LaTex','FontSize', 20)
legend('j = 82', 'j = 163')

figure
p(3) = imagesc(standardized_cost_matrix);
c = colorbar;
ylabel(c, '$\hat{C_i}(\bar{x_j})$', 'Interpreter', 'LaTex','FontSize', 20)
title('Imagesc Plot Of Normalized Cost', 'Interpreter', 'LaTex','FontSize', 20)
xlabel('j = Modes', 'Interpreter', 'LaTex','FontSize', 20)
ylabel('i = Experiments', 'Interpreter', 'LaTex','FontSize', 20)

[mem, father, costs] = KBestModes(standardized_cost_matrix, 243);
min_cost_per_k_modes = zeros(243, 1);
for  i = 1:243
    min_cost_per_k_modes(i) = min(costs(:, i));
end
figure
p(4) = plot(min_cost_per_k_modes);
title('Plot Of Normalized Total Cost Over Number Of Modes', 'Interpreter', 'LaTex','FontSize', 20)
xlabel('Number Of Modes In The Family', 'Interpreter', 'Latex','FontSize', 20)
ylabel('$ \bar{C}(\bar{x}; m_1, ..., m_K)$', 'Interpreter', 'LaTex','FontSize', 20)
