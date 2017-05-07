%% Clear
clear all;
close all;
clc;
%% Setup
run('Setup.m');
import Models.QSPusherSlider
import MPCSolvers.FOMSolver
import Combinatorics.VariationWithRepetition
DateString = datestr(clock, 'mm_dd_yyyy__HH_MM_SS');
%% Optimization Hybrid States Set Up
c2 = QSPusherSlider.c^2;
hybrid_states_map = horzcat(PusherSliderStates.QSSticking(QSPusherSlider.a, QSPusherSlider.nu_pusher, c2), ...
                            PusherSliderStates.QSSlidingUp(QSPusherSlider.a, QSPusherSlider.nu_pusher, c2), ...
                            PusherSliderStates.QSSlidingDown(QSPusherSlider.a, QSPusherSlider.nu_pusher, c2));
c2_pert = QSPusherSlider.c_pert^2;
%% Euler Integration Hybrid States Set Up
real_states_map = horzcat(PusherSliderStates.QSSticking(QSPusherSlider.a, QSPusherSlider.nu_pusher_pert, c2_pert), ...
                          PusherSliderStates.QSSlidingUp(QSPusherSlider.a, QSPusherSlider.nu_pusher_pert, c2_pert), ...
                          PusherSliderStates.QSSlidingDown(QSPusherSlider.a, QSPusherSlider.nu_pusher_pert, c2_pert));
%% Hybrid Modes Parameters and Setup
fom_steps = 5;
% ext_fom_steps = 35;
% diff_fom_steps = ext_fom_steps - fom_steps;
hybrid_modes = VariationWithRepetition(length(hybrid_states_map), fom_steps);
% hybrid_modes = [hybrid_modes, ones(size(hybrid_modes, 1), diff_fom_steps)];
%% MPC Parameters and Setup
Q_MPC = 1 * diag([50,50,.1,0]); % 10 * diag([1,10,10,0]);
Q_MPC_final = 100 * diag([1,1,.1,0]); % 200 * diag([.1,10,1,0]);
R_MPC = 0.5 * diag([1,1]);
u_lower_bound = [-0.01; 0.1];
u_upper_bound = [0.1; 0.1];
x_lower_bound = [1; 1; 20; QSPusherSlider.a];
x_upper_bound = [1; 1; 20; QSPusherSlider.a];
h_opt = 0.03;
solver = FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes, 0); % FOM without chameleon
%% Solve MPC and Integrate
tf = fom_steps * h_opt; % We only take one step of the Euler integration
u_thrust = 0.05;
u_star = @(t)([u_thrust * ones(size(t)); 0 .*t]);
x_star_line = @(t)([u_thrust .* t; 0.*t; 0.*t; 0.*t]); % For the straight line
c = QSPusherSlider.c;
px = QSPusherSlider.a / 2.0;
A = @(d)(u_thrust / (c^2 + px^2 + d^2));
x_star_circle = @(d)(@(t)([-(sin(-(d*A(d)).*t) * (c^2 + px^2) + cos(-(d*A(d)) .* t) * px*d - px*d) / d; ...
    (cos(-(d*A(d)).*t) * (c^2 + px^2) - sin(-(d*A(d)).* t) * px*d - c^2 - px^2) / d; ...
    -(d*A(d)) .* t; ... 
    d * ones(size(t))])); % For the circle

number_of_curvature_pairs = 10;
d = linspace(-3 * QSPusherSlider.b/4.0, 3 * QSPusherSlider.b/4.0, 2 * number_of_curvature_pairs);
eps = 0.15; % TODO: Automatically generate it
angular_eps = pi / 2; % 90% in each direction. Otherwise we choose another side
sim_number = 1000; % Number of simulations
n_rand = @(n, m, eps)((rand(n, m) - 0.5) * eps);
x0 = [n_rand(sim_number, 2, 2 * eps), n_rand(sim_number, 1, 2 * angular_eps), n_rand(sim_number, 1, QSPusherSlider.a)].';
path_name = 'SimulationResults/GeneralTrajectoryLearning';
if  exist(path_name, 'dir') ~= 7
    mkdir(pwd, path_name);
end

normalize_matrix = @(A)((A - nanmean(A, 2) * ones(1, size(A, 2))) ./ (nanstd(A, 0, 2) * ones(1, size(A, 2))));
% number_of_modes = 20;
best_modes = cell(1, 2 * number_of_curvature_pairs + 1);
best_costs = cell(1, 2 * number_of_curvature_pairs + 1);
for k = 0:2 * number_of_curvature_pairs
% for k = 14
    disp(k)
    if k == 0
        x_star = x_star_line;
    else
        x_star = x_star_circle(d(k));
    end
    costs_matrix = zeros(sim_number, size(hybrid_modes, 1));
    for i = 1:sim_number
        disp(i)
        for j = 1:size(hybrid_modes, 1)
            op = solver.GetOptimizationProblem(0, x_star, u_star, x0(:,i), hybrid_modes(j, :));
            try
                [~, ~, costs_matrix(i, j)] = op.solve;
            catch
                disp('Infeasible');
                costs_matrix(i, j) = nan;
            end
        end
    end
    mkdir(pwd, strcat(path_name, '/', DateString));
    if k == 0
        save(strcat(path_name, '/', DateString, '/', 'cost_data_for_simnumber_', int2str(sim_number), '_straight_line'), 'costs_matrix');
    else
        save(strcat(path_name, '/', DateString, '/', 'cost_data_for_simnumber_', int2str(sim_number), '_curvature_dx1000_', int2str(floor(d(k) * 1000))), 'costs_matrix');
    end
    standardized_cost_matrix = normalize_matrix(costs_matrix);
    standardized_cost_matrix(isnan(standardized_cost_matrix)) = Inf;
    if k == 0
        save(strcat(path_name, '/', DateString, '/', 'std_cost_data_for_simnumber_', int2str(sim_number), '_straight_line'), 'standardized_cost_matrix');
    else
        save(strcat(path_name, '/', DateString, '/', 'std_cost_data_for_simnumber_', int2str(sim_number), '_curvature_dx1000_', int2str(floor(d(k) * 1000))), 'standardized_cost_matrix');
    end
    [best_costs{k+1}, best_modes{k+1}] = sort(nanmean(standardized_cost_matrix, 1), 'ascend');
end