import Models.QSPusherSlider
folder = '04_28_2017__11_59_32';
best_modes = zeros(10, 21);
path_name = 'SimulationResults/GetBestModes';
if  exist(path_name, 'dir') ~= 7
    mkdir(pwd, path_name);
end
DateString = datestr(clock, 'mm_dd_yyyy__HH_MM_SS');
number_of_curvature_pairs = 10;
d = linspace(-3 * QSPusherSlider.b/4.0, 3 * QSPusherSlider.b/4.0, 2 * number_of_curvature_pairs);
int_d = int2str(floor(d * 1000));

% figure
hold on
for k = 1:21
    if k == 21
        load(strcat('SimulationResults/GeneralTrajectoryLearning/', folder, '/cost_data_for_simnumber_1000_straight_line.mat'))
        load(strcat('SimulationResults/GeneralTrajectoryLearning/', folder, '/std_cost_data_for_simnumber_1000_straight_line.mat'))
    else
        load(strcat('SimulationResults/GeneralTrajectoryLearning/', folder, '/cost_data_for_simnumber_1000_curvature_dx1000_', int2str(floor(d(k) * 1000)), '.mat'))
        load(strcat('SimulationResults/GeneralTrajectoryLearning/', folder, '/std_cost_data_for_simnumber_1000_curvature_dx1000_', int2str(floor(d(k) * 1000)), '.mat'))
    end
    index = zeros(243, 1);
    cost_greedy = inf * ones(243, 1);
    vec = cell(243, 1);
    cost_greedy(1) = mean(standardized_cost_matrix(:,1));
    index(1) = 1;
    vec{1} = standardized_cost_matrix(:, index(1));
    for i = 2:243
        for j = 1:243
            loc_vec = min(vec{i-1}, standardized_cost_matrix(:, j));
            if mean(loc_vec) < cost_greedy(i)
                cost_greedy(i) = mean(loc_vec);
                index(i) = j;
                vec{i} = loc_vec;
            end
        end
    end
    best_modes(:, k) = index(1:10);
%     plot(cost_greedy,'LineWidth',2);
end
% leg = {'Straigth Line'};
% leg = [leg int2str((1:20).')];
legend()
mkdir(pwd, strcat(path_name, '/', DateString));
save(strcat(path_name, '/', DateString, '/best_modes'), 'best_modes', 'd');