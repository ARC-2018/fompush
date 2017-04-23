%% Clear
clear all;
close all;
clc;
%% Setup
run('Setup.m');
import Models.QSPusherSlider
import MPCSolvers.FOMSolver
import Combinatorics.VariationWithRepetition
import LPusher.Simulator2Controller
import LPusher.Controller2Simulator
DateString = datestr(clock, 'mm_dd_yyyy__HH_MM_SS');
%% Optimization Hybrid States Set Up
c2 = QSPusherSlider.c^2;
hybrid_states_map = horzcat(LPusher.LPStickingCont(), ...
                            LPusher.LPUpCont(), ...
                            LPusher.LPDownCont());
c2_pert = QSPusherSlider.c_pert^2;
%% Euler Integration Hybrid States Set Up
real_states_map = LPusher.LPSimulator;
%% Hybrid Modes Parameters and Setup
fom_steps = 5;
ext_fom_steps = 35;
diff_fom_steps = ext_fom_steps - fom_steps;
all_modes = VariationWithRepetition(length(hybrid_states_map), fom_steps);
% 100 193 198
hybrid_modes1 = [all_modes(1,:), ones(1, diff_fom_steps);all_modes(82,:), ones(1, diff_fom_steps);all_modes(163,:), ones(1, diff_fom_steps)];
hybrid_modes2 = [all_modes(100,:), ones(1, diff_fom_steps);all_modes(193,:), ones(1, diff_fom_steps);all_modes(198,:), ones(1, diff_fom_steps)];
%% MPC Parameters and Setup
[~, Q_MPC_final, ~, ~, ~, ~, ~, h_opt, h_step] = GetTestParams();
Q_MPC = 10*diag([1,3,1,0]);
R_MPC = 1*diag([1,1,1,1,0]);
x_lower_bound = [1000;1000;1000;1000];
x_upper_bound = [1000;1000;1000;1000];
u_lower_bound = [-0.01;-0.01;1000;1000;1000];
u_upper_bound = [1000;1000;1000;1000;1000];
solvers{1} = FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, @Simulator2Controller, @Controller2Simulator, hybrid_modes1, 1); % FOM with chameleon
solvers{2} = FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, @Simulator2Controller, @Controller2Simulator, hybrid_modes2, 1);
%% Solve MPC and Integrate
tf = 1; % We only take one step of the Euler integration
u_thrust = 0.05;
x_s = @(t)([u_thrust .* t 0.* t 0.* t u_thrust .* t - QSPusherSlider.a/2.0 0.* t 0.* t]');
u_s = @(t)([u_thrust * ones(size(t)) 0.* t 0.* t]');
%% Initial conditions
x_c_0 = [0.;0.05;15*pi/180;0]; 
[x_s_0, ~] = LPusher.Controller2Simulator(x_c_0);
euler_integrator = EulerIntegration(real_states_map, h_step);

cost = cell(length(solvers), 1);
for i = 1:length(solvers)
    animator = LPusher.Animator('test');
    animator.NumSim = 1;
    tic;
    [x_state, u_state, x_bar, u_bar, t, modes, costs, animator] = euler_integrator.IntegrateMPC(0, tf, x_s_0, x_s, u_s, solvers{i}, animator);
    toc;
    x_star_realized = x_state - x_bar;
    animator.plot(0);
    animator.plot(1)
    animator.plot(2)
    animator.Animate(1)
%     cost{i} = dot(x_bar, Q_MPC * x_bar) + dot(u_bar, R_MPC * u_bar);
end
% for i = 1:length(solvers) 
%     sum(cost{i})
% end


