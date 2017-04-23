%% Setup
run('Setup.m');
import MPCSolvers.FOMSolver
import Combinatorics.VariationWithRepetition
import LPusher.LPusherSlider
import LPusher.Friction
DateString = datestr(clock, 'mm_dd_yyyy__HH_MM_SS');
%% Optimization Hybrid States Set Up
hybrid_states_map = horzcat(LPusher.LPStickingCont(), ...
                            LPusher.LPUpCont(), ...
                            LPusher.LPDownCont());
%% Euler Integration Hybrid States Set Up
real_states_map = LPusher.LPSimulator;
%% Hybrid Modes Parameters and Setup
fom_steps = 5;
hybrid_modes = VariationWithRepetition(length(hybrid_states_map), fom_steps);
% ext_fom_steps = 35; % FOR LEARNING WITH TAIL
% diff_fom_steps = ext_fom_steps - fom_steps;
% hybrid_modes = [hybrid_modes ones(size(hybrid_modes, 1), diff_fom_steps)];
%% MPC Parameters and Setup
Q_MPC_final = 200 * diag([1,1,.1,0]); % 200 * diag([.1,10,1,0]);
Q_MPC = 10*diag([1,3,1,0]);
R_MPC = 1*diag([1,1,1,1,0]);
x_lower_bound = [1000;1000;1000;1000];
x_upper_bound = [1000;1000;1000;1000];
u_lower_bound = [-0.01;-0.01;1000;1000;1000];
u_upper_bound = [1000;1000;1000;1000;1000];
solver = FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, @Simulator2Controller, @Controller2Simulator, hybrid_modes, 0); % FOM without chameleon
%% Solve MPC and Integrate
tf = fom_steps * h_opt; % We only take one step of the Euler integration
u_thrust = 0.05;
u_s_star_line = @(t)([u_thrust * ones(size(t)) 0.* t 0.* t]');
x_s_star_line = @(t)([u_thrust .* t 0.* t 0.* t u_thrust .* t - LPusherSlider.a/2.0 0.* t 0.* t]');% For the straight line

% u_thrust = 0.05;
% x_s_star_line = @(t)([u_thrust .* t 0.* t 0.* t u_thrust .* t - LPusherSlider.a/2.0 0.* t 0.* t]');
% u_s_star_line = @(t)([u_thrust * ones(size(t)) 0.* t 0.* t]');
% c = LPusherSlider.c;
% px = LPusherSlider.a / 2.0;
% A = @(d)(u_thrust / (c^2 + px^2 + d^2));
% x_s_star_circle = @(d)(@(t)([-(sin(-(d*A(d)).*t) * (c^2 + px^2) + cos(-(d*A(d)) .* t) * px*d - px*d) / d; ...
%     (cos(-(d*A(d)).*t) * (c^2 + px^2) - sin(-(d*A(d)).* t) * px*d - c^2 - px^2) / d; ...
%     -(d*A(d)) .* t; ... 
%     d * ones(size(t))])); % For the circle

number_of_curvature_pairs = 0;
d = linspace(-3 * LPusherSlider.b/4.0, 3 * LPusherSlider.b/4.0, 2 * number_of_curvature_pairs);
eps = 0.5; % TODO: Automatically generate it
angular_eps = pi / 4; % 45% in each direction. Otherwise we choose another side
sim_number = 1000; % Number of simulations
n_rand = @(n, m, eps)((rand(n, m) - 0.5) * eps);

x_c_0 = [n_rand(sim_number, 2, 2 * eps), n_rand(sim_number, 1, 2 * angular_eps), n_rand(sim_number, 1, LPusherSlider.a - LPusherSlider.lp)].';
x_s_0 = zeros(6, sim_number);
for i = 1:size(x_c_0, 2) %TODO: MAKE FASTER
    x_s_0(:, i) = LPusher.Controller2Simulator(x_c_0(:, i));
end

path_name = 'SimulationResults/LinePusherLearning';
if  exist(path_name, 'dir') ~= 7
    mkdir(pwd, path_name);
end

normalize_matrix = @(A)((A - nanmean(A, 2) * ones(1, size(A, 2))) ./ (nanstd(A, 0, 2) * ones(1, size(A, 2))));
% number_of_modes = 20;
best_modes = cell(1, 2 * number_of_curvature_pairs + 1);
best_costs = cell(1, 2 * number_of_curvature_pairs + 1);
for k = 0:2 * number_of_curvature_pairs
    disp(k)
    if k == 0
        x_s_star = x_s_star_line;
        u_s_star = u_s_star_line;
    else
        x_s_star = x_s_star_circle(d(k));
        u_s_star = u_s_star_circle(d(k));
    end
    costs_matrix = zeros(sim_number, size(hybrid_modes, 1));
    for i = 1:sim_number
        disp(i)
        for j = 1:size(hybrid_modes, 1)
            op = solver.GetOptimizationProblem(0, x_s_star, u_s_star, x_s_0(:,i), hybrid_modes(j, :));
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