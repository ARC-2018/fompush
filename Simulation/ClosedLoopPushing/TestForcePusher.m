%% Clear
% clear all;
close all;
clc;
%% Setup
run('Setup.m');
import Models.QSPusherSlider
import MPCSolvers.Force1pointSolver
import Combinatorics.VariationWithRepetition
%% Hybrid Modes Parameters and Setup
fom_steps = 5;
ext_fom_steps = 35;
diff_fom_steps = ext_fom_steps - fom_steps;
all_modes = VariationWithRepetition(3, fom_steps);
% h_test = [1 1; 2 1; 3 1];
hybrid_modes1 = [all_modes(1,:), ones(1, diff_fom_steps);all_modes(82,:), ones(1, diff_fom_steps);all_modes(163,:), ones(1, diff_fom_steps)];
hybrid_modes2 = [all_modes(1,:), ones(1, diff_fom_steps);all_modes(203,:), ones(1, diff_fom_steps);all_modes(160,:), ones(1, diff_fom_steps)];
hybrid_modes3 = [hybrid_modes2; all_modes(82,:), ones(1, diff_fom_steps)];
hybrid_modes4 = [hybrid_modes3; all_modes(163,:), ones(1, diff_fom_steps);];
hybrid_modes5 = [hybrid_modes4;all_modes(235,:), ones(1, diff_fom_steps)];
hybrid_modes6 = [hybrid_modes5;all_modes(235,:), ones(1, diff_fom_steps)];
hybrid_modes7 = [hybrid_modes6;all_modes(162,:), ones(1, diff_fom_steps)];


Q_MPC = 10*diag([1,3,1,0]);
Q_MPC_final = 2000*diag([1,3,1,0]);
R_MPC = 1*diag([1,1,0]);

% Q_MPC = 1 * diag([50,50,.1,0]);
% Q_MPC_final = 100 * diag([1,1,.1,0]); % 200 * diag([.1,10,1,0]);
% R_MPC = 0.5 * diag([1,1,1]);
u_lower_bound = [100; 100; 100]; % TODO: Change
u_upper_bound = [100; 100; 100]; % TODO: Change
x_lower_bound = [1; 1; 20; QSPusherSlider.a/2];
x_upper_bound = [1; 1; 20; QSPusherSlider.a/2];
h_opt = 0.03;
h_step = 0.01;
solvers = {};
% solvers = [solvers, {Force1pointSolver(Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, h_test)}];
solvers = [solvers, {Force1pointSolver(Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes2, 1)}];
% solvers = [solvers,{Force1pointSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes3, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes4, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes5, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes6, 1)}];
% solvers = [solvers,{FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes7, 1)}];
euler_integrator = ForceEulerIntegration(Oneptpusher.Simulator(), h_step);
%% Solve MPC and Integrate
t0 = 0;
tf = 4;
u_thrust = 0.05;
uns = spline([1 2], [u_thrust u_thrust; 0 0]);
xns = spline([0 1], [0 u_thrust; 0 0; 0 0; -QSPusherSlider.a/2 u_thrust - QSPusherSlider.a/2; 0 0]);

xsStarTemp = @(t)([u_thrust*obj.t_star(lv1) 0 0 u_thrust*obj.t_star(lv1)-obj.a/2 0]');
usStarTemp = @(t)([u_thrust 0]');

% x_star = @(t)([u_thrust .* t; 0.*t; 0.*t; 0.*t]); % For the straight line
% d = 0.5 * QSPusherSlider.b/2.0
% d = 0.0249;
% c = QSPusherSlider.c;
% px = QSPusherSlider.a / 2.0;
% A = u_thrust / (c^2 + px^2 + d^2);
% x_star = @(t)([-(sin(-(d*A).*t) * (c^2 + px^2) + cos(-(d*A) .* t) * px*d - px*d) / d; (cos(-(d*A).*t) * (c^2 + px^2) - sin(-(d*A).* t) * px*d - c^2 - px^2) / d; -(d*A) .* t; d * ones(size(t))]); % For the circle

% xc_0{1} = [0.;0.05;0*pi/180;0]; 
% xs_0{1} = c.coordinateTransformCS(xc_0{1});
xs0 =  [0; 0.05; 0 * pi / 180; -QSPusherSlider.a / 2; 0.05];
% xs0 = [0, .1, 30 * pi / 180, 0, -QSPusherSlider.a/2, .1]; % For the straight line
% x0 = [0.0, -0.05, -0.1, d]; % For the circle

% eps = 0.5; % TODO: Automatically generate it
% angular_eps = pi / 2; % 90% in each direction. Otherwise we choose another side
% sim_number = 5; % Number of simulations
% n_rand = @(n, m, eps)((rand(n, m) - 0.5) * eps);
% x0 = [n_rand(sim_number, 2, 2 * eps), n_rand(sim_number, 1, 2 * angular_eps), n_rand(sim_number, 1, QSPusherSlider.a)];
% x0 = [0.03, .05, -0.2, d]; % high circle disturbance USE IN THESIS
%% Animation Parameters
animator = Animation.Animator(QSPusherSlider());
path_name = 'SimulationResults/TestForcePusher';
if  exist(path_name, 'dir') ~= 7
    mkdir(pwd, path_name);
end
video_names = {strcat(path_name, '/fom_3_modes'), strcat(path_name, '/fom_4_modes'), ...
    strcat(path_name, '/fom_5_modes'), strcat(path_name, '/fom_6_modes'), strcat(path_name, '/fom_7_modes'), strcat(path_name, '/fom_8_modes'), strcat(path_name, '/fom_9_modes')};
costs = cell(length(solvers), 1);
xs_bars = cell(length(solvers), 1);
xs = cell(length(solvers), 1);
xc = cell(length(solvers), 1);
us_bars = cell(length(solvers), 1);
fh_chosen_x = cell(length(solvers), 1);
fh_considered_x = cell(length(solvers), 1);
f_horizon_u_bar = cell(length(solvers), 1);
for i = 1:length(solvers)
%     for j = 1:sim_number
    [xs{i}, xc{i}, us, uc, xs_bars{i}, us_bars{i}, t, costs{i}, fh_chosen_x{i}, fh_considered_x{i}] = euler_integrator.IntegrateMPC(t0, tf, xs0, xns, uns, solvers{i});
%     end
%     global_cost{i} = dot(x_bars{i}, Q_MPC * x_bars{i}) + dot(u_bars{i}, R_MPC * u_bars{i});
end
% for i = 1:length(solvers) 
%     sum(global_cost{i})
% end
for i = 1:length(solvers)
    video_name = video_names{i};
    frame_rate = length(xs{i}(1, :)) / ((tf - t0));
    animator.AnimateTrajectory(frame_rate, video_name, solvers{i}.GetXC(xns, t.'), xc{i}, fh_chosen_x{i}, fh_considered_x{i}, 0);
end
% h = figure;
% hold on;
% colors = {'blue',[1 .5 0],'green','black', 'yellow', 'red'};
% lstile = {':','-', '--',':','-','-','--'};
% for i = 1:length(solvers)
%     cost = zeros(1, size(x_bars{i}, 2));
%     for j = 1:size(x_bars{i}, 2)
%         cost(j) = x_bars{i}(:, j).' * Q_MPC * x_bars{i}(:, j) + u_bars{i}(:, j).' * R_MPC * u_bars{i}(:, j);
%     end
%     plot(t, cost, 'color', colors{i}, 'LineStyle', lstile{i}, 'LineWidth', 2);
% end
% l = legend('3 original modes', '3 best modes', '4 best modes', '5 best modes', '7 best modes', 'best 8 modes', 'Interpreter', 'Latex');
% set(l,'FontSize',20);
% title('Local cost over time', 'Interpreter', 'LaTex','FontSize', 30)
% xlabel('Time', 'Interpreter', 'LaTex','FontSize', 30)
% ylabel('Value of $\bar{x}_i^T Q \bar{x}_i + \bar{u}_i^T R \bar{u}_i$', 'Interpreter', 'LaTex','FontSize', 30)
% h = figure;
% hold on;
% % for i = [1:length(solvers)]
% for i = [1:length(solvers)]
%     costs{i};
%     plot(t, costs{i},'color',colors{i}, 'LineStyle', lstile{i},'LineWidth',2);
% end
% % legend('3 modes', '3 best modes', '7 best modes');
% l2 = legend('original 3 modes', '3 best modes', '4 best modes', '5 best modes', '6 best modes', '7 best modes');
% set(l2,'FontSize',20);
% title('Local MPC cost over time', 'Interpreter', 'LaTex','FontSize', 30)
% xlabel('Time', 'Interpreter', 'LaTex','FontSize', 30)
% ylabel('Value of $C_i$', 'Interpreter', 'LaTex','FontSize', 30)
