%% Clear
clear all;
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
ext_fom_steps = 35;
diff_fom_steps = ext_fom_steps - fom_steps;
all_modes = VariationWithRepetition(length(hybrid_states_map), fom_steps);
hybrid_modes1 = [all_modes(1,:), ones(1, diff_fom_steps);all_modes(82,:), ones(1, diff_fom_steps);all_modes(163,:), ones(1, diff_fom_steps)];
hybrid_modes2 = [hybrid_modes1;all_modes(28,:), ones(1, diff_fom_steps)];
hybrid_modes3 = [hybrid_modes2;all_modes(10,:), ones(1, diff_fom_steps)];


hybrid_modes4 = [hybrid_modes3;all_modes(4,:), ones(1, diff_fom_steps)];
hybrid_modes5 = [hybrid_modes4;all_modes(2,:), ones(1, diff_fom_steps)];
hybrid_modes6 = [hybrid_modes5;all_modes(3,:), ones(1, diff_fom_steps)];
hybrid_modes7 = [hybrid_modes6;all_modes(7,:), ones(1, diff_fom_steps)];
[Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, h_step] = GetTestParams();
Q_MPC = 1 * diag([50,50,.1,0]);

solvers = {FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes1, 1)}; % FOM with chameleon
solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes2, 1)}];
solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes3, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes4, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes5, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes6, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes7, 1)}];
euler_integrator = EulerIntegration(real_states_map, h_step);
%% Solve MPC and Integrate
t0 = 0;
tf = 10;
u_thrust = 0.05;
u_star = @(t)([u_thrust * ones(size(t)); 0 .*t]); % To ensure it can be substituted by a vector
% x_star = @(t)([u_thrust .* t; 0.*t; 0.*t; 0.*t]); % For the straight line
d = 0.5 * QSPusherSlider.b/2.0;
c = QSPusherSlider.c;
px = QSPusherSlider.a / 2.0;
A = u_thrust / (c^2 + px^2 + d^2);
x_star = @(t)([-(sin(-(d*A).*t) * (c^2 + px^2) + cos(-(d*A) .* t) * px*d - px*d) / d; (cos(-(d*A).*t) * (c^2 + px^2) - sin(-(d*A).* t) * px*d - c^2 - px^2) / d; -(d*A) .* t; d * ones(size(t))]); % For the circle
% x0 = [0, .05, 30 * pi / 180, 0]; % For the straight line
x0 = [0.0, -0.03, -0.1, d]; % For the circle
% x0 = [0.03, .15, -0.2, d]; % high circle disturbance USE IN THESIS
%% Animation Parameters
animator = Animation.Animator(QSPusherSlider());
path_name = 'SimulationResults/TestLineWithOptimalModes';
if  exist(path_name, 'dir') ~= 7
    mkdir(pwd, path_name);
end
video_names = {strcat(path_name, '/fom_3_modes'), strcat(path_name, '/fom_4_modes'), ...
    strcat(path_name, '/fom_5_modes'), strcat(path_name, '/fom_6_modes'), strcat(path_name, '/fom_7_modes'), strcat(path_name, '/fom_8_modes'), strcat(path_name, '/fom_9_modes')};
cost = cell(length(solvers), 1);
x_bars = cell(length(solvers), 1);
x_states = cell(length(solvers), 1);
u_bars = cell(length(solvers), 1);
for i = 1:length(solvers)
    tic;
    [x_states{i}, u_state, x_bars{i}, u_bars{i}, t, modes, costs] = euler_integrator.IntegrateMPC(t0, tf, x0, x_star, u_star, solvers{i});
    toc;
    cost{i} = dot(x_bars{i}, Q_MPC * x_bars{i}) + dot(u_bars{i}, R_MPC * u_bars{i});
end
for i = 1:length(solvers) 
    sum(cost{i})
end
for i = 1:length(solvers)
    video_name = video_names{i};
    frame_rate = length(x_states{i}(1, :)) / (50 * (tf - t0));
    animator.AnimateTrajectory(frame_rate, video_name, x_states{i} - x_bars{i}, x_states{i});
end
h = figure;
hold on;
for i = [1, 2, 3]%, 4, 7]
    cost = zeros(1, size(x_bars{i}, 2));
    for j = 1:size(x_bars{i}, 2)
        cost(j) = x_bars{i}(:, j).' * Q_MPC * x_bars{i}(:, j) + u_bars{i}(:, j).' * R_MPC * u_bars{i}(:, j);
    end
    plot(t, cost);
end
legend('3 modes', '4 modes', '5 modes', '9 modes');
title('Local cost over time', 'Interpreter', 'LaTex','FontSize', 20)
xlabel('Time', 'Interpreter', 'LaTex','FontSize', 20)
ylabel('Value of $C_i$', 'Interpreter', 'LaTex','FontSize', 20)
