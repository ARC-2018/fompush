%% Clear
% clear all;
close all;
clc;
%% Setup
run('Setup.m');
import Models.QSPusherSlider
import MPCSolvers.FOMSolver
import MPCSolvers.MIQPSolver
import Combinatorics.VariationWithRepetition
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
miqp_steps = 5;
fom_steps = 5;
clustering_factor = 7;
ext_fom_steps = 35;
diff_fom_steps = ext_fom_steps - fom_steps;
all_modes = VariationWithRepetition(length(hybrid_states_map), fom_steps);
hybrid_modes1 = [all_modes(1,:), ones(1, diff_fom_steps);all_modes(82,:), ones(1, diff_fom_steps);all_modes(163,:), ones(1, diff_fom_steps)];
hybrid_modes2 = [all_modes(1,:), ones(1, diff_fom_steps);all_modes(203,:), ones(1, diff_fom_steps);all_modes(160,:), ones(1, diff_fom_steps)];
hybrid_modes3 = [hybrid_modes2; all_modes(82,:), ones(1, diff_fom_steps)];
hybrid_modes4 = [hybrid_modes3; all_modes(163,:), ones(1, diff_fom_steps);];
hybrid_modes5 = [hybrid_modes4;all_modes(235,:), ones(1, diff_fom_steps)];
hybrid_modes6 = [hybrid_modes5;all_modes(235,:), ones(1, diff_fom_steps)];
hybrid_modes7 = [hybrid_modes6;all_modes(162,:), ones(1, diff_fom_steps)];
[Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, h_step] = GetTestParams();
Q_MPC = 1 * diag([50,50,.1,0]);
solvers = {}
% solvers = {FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, [1, 1, 1, ones(1, 10)], 0)};
% solvers = [solvers, {FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, [1, 1, ones(1, 10),;3, 2, ones(1, 10)], 1)}];
% solvers = [solvers, {FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, [1, 1, ones(1, 10),;3, 2, ones(1, 10); 2, 1, ones(1, 10)], 1)}];
% solvers = [solvers, {MIQPSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, miqp_steps, 1)}]; % MIQP
% solvers = [solvers, {MIQPSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, miqp_steps * clustering_factor, clustering_factor)}]; % MIQP with clustering
% solvers = [solvers, {FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes1, 1)}]; % FOM with chameleon
solvers = [solvers, {FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes1, 1)}]; % FOM with chameleon
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes2, 1)}];
solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes3, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes4, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes5, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes6, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes7, 1)}];
euler_integrator = EulerIntegration(real_states_map, h_step);
%% Solve MPC and Integrate
t0 = 0;
tf = 15;
u_thrust = 0.05;
u_star = @(t)([u_thrust * ones(size(t)); 0 .*t]); % To ensure it can be substituted by a vector
% x_star = @(t)([u_thrust .* t; 0.*t; 0.*t; 0.*t]); % For the straight line
% d = 0.5 * QSPusherSlider.b/2.0
d = 0.0249;
c = QSPusherSlider.c;
px = QSPusherSlider.a / 2.0;
A = u_thrust / (c^2 + px^2 + d^2);
x_star = @(t)([-(sin(-(d*A).*t) * (c^2 + px^2) + cos(-(d*A) .* t) * px*d - px*d) / d; (cos(-(d*A).*t) * (c^2 + px^2) - sin(-(d*A).* t) * px*d - c^2 - px^2) / d; -(d*A) .* t; d * ones(size(t))]); % For the circle
% x0 = [0, .1, 30 * pi / 180, 0]; % For the straight line
x0 = [0.0, -0.05, -0.1, d]; % For the circle

% eps = 0.5; % TODO: Automatically generate it
% angular_eps = pi / 2; % 90% in each direction. Otherwise we choose another side
% sim_number = 5; % Number of simulations
% n_rand = @(n, m, eps)((rand(n, m) - 0.5) * eps);
% x0 = [n_rand(sim_number, 2, 2 * eps), n_rand(sim_number, 1, 2 * angular_eps), n_rand(sim_number, 1, QSPusherSlider.a)];
% x0 = [0.03, .05, -0.2, d]; % high circle disturbance USE IN THESIS
%% Animation Parameters
animator = Animation.Animator(QSPusherSlider());
path_name = 'SimulationResults/TestCurveWithDPLearning';
if  exist(path_name, 'dir') ~= 7
    mkdir(pwd, path_name);
end
video_names = {strcat(path_name, '/fom_3_modes'), strcat(path_name, '/fom_4_modes'), ...
    strcat(path_name, '/fom_5_modes'), strcat(path_name, '/fom_6_modes'), strcat(path_name, '/fom_7_modes'), strcat(path_name, '/fom_8_modes'), strcat(path_name, '/fom_9_modes')};
costs = cell(length(solvers), 1);
x_bars = cell(length(solvers), 1);
x_states = cell(length(solvers), 1);
u_bars = cell(length(solvers), 1);
fh_chosen_x = cell(length(solvers), 1);
fh_considered_x = cell(length(solvers), 1);
f_horizon_u_bar = cell(length(solvers), 1);
for i = 1:length(solvers)
%     for j = 1:sim_number
    [x_states{i}, u_state, x_bars{i}, u_bars{i}, t, modes, costs{i}, fh_chosen_x{i}, fh_considered_x{i}] = euler_integrator.IntegrateMPC(t0, tf, x0, x_star, u_star, solvers{i});
%     end
    global_cost{i} = dot(x_bars{i}, Q_MPC * x_bars{i}) + dot(u_bars{i}, R_MPC * u_bars{i});
end
for i = 1:length(solvers) 
    sum(global_cost{i})
end
for i = 1:length(solvers)
    video_name = video_names{i};
    frame_rate = length(x_states{i}(1, :)) / ((tf - t0));
    animator.AnimateTrajectory(frame_rate, video_name, x_states{i} - x_bars{i}, x_states{i}, fh_chosen_x{i}, fh_considered_x{i}, 0);
end
h = figure;
hold on;
colors = {'blue',[1 .5 0],'green','black', 'yellow', 'red'};
lstile = {':','-', '--',':','-','-','--'};
for i = 1:length(solvers)
    cost = zeros(1, size(x_bars{i}, 2));
    for j = 1:size(x_bars{i}, 2)
        cost(j) = x_bars{i}(:, j).' * Q_MPC * x_bars{i}(:, j) + u_bars{i}(:, j).' * R_MPC * u_bars{i}(:, j);
    end
    plot(t, cost, 'color', colors{i}, 'LineStyle', lstile{i}, 'LineWidth', 2);
end
l = legend('3 original modes', '3 best modes', '4 best modes', '5 best modes', '7 best modes', 'best 8 modes', 'Interpreter', 'Latex');
set(l,'FontSize',20);
title('Local cost over time', 'Interpreter', 'LaTex','FontSize', 30)
xlabel('Time', 'Interpreter', 'LaTex','FontSize', 30)
ylabel('Value of $\bar{x}_i^T Q \bar{x}_i + \bar{u}_i^T R \bar{u}_i$', 'Interpreter', 'LaTex','FontSize', 30)
h = figure;
hold on;
% for i = [1:length(solvers)]
for i = [1:length(solvers)]
    costs{i};
    plot(t, costs{i},'color',colors{i}, 'LineStyle', lstile{i},'LineWidth',2);
end
% legend('3 modes', '3 best modes', '7 best modes');
l2 = legend('original 3 modes', '3 best modes', '4 best modes', '5 best modes', '6 best modes', '7 best modes');
set(l2,'FontSize',20);
title('Local MPC cost over time', 'Interpreter', 'LaTex','FontSize', 30)
xlabel('Time', 'Interpreter', 'LaTex','FontSize', 30)
ylabel('Value of $C_i$', 'Interpreter', 'LaTex','FontSize', 30)
