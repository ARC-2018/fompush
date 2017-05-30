classdef ForceEulerIntegration
    % Class to implement Euler Integration for certain desired
    % trajectories. It currently implements simple integration and MPC
    % solving
properties (Access = private)
    h_step;
    simulator;
end
methods
function obj = ForceEulerIntegration(simulator, h_step)
    obj.h_step = h_step;
    obj.simulator = simulator;
end

function [] = GetScoresOneStep(obj, t0, x_star, u_star, mpc_solver) % TODO
    local_mpc_solver = mpc_solver; % To avoid corrupting the mpc starting data
    assert(length(obj.real_states_map(1).x) == length(x0), 'x_state in Euler integration and real states has different size');
    x_state = zeros(length(x0), 1);
    x_state(:, 1) = x0; 
    number_of_controllers = length(u_star(t0));
    u_bar = zeros(number_of_controllers, 1);
        [modes, costs] = local_mpc_solver.SolveMPC(t0, x_star, u_star, x_state(:, step)); % Solve Controller MPC problem and get control input
        for real_state_index = 1:(length(obj.real_states_map) + 1) % Determine which state this controller actually falls into
            assert(real_state_index < length(obj.real_states_map) + 1, sprintf('At step %d of EulerIntegration, the real state to apply Euler integration cannot be determined', step));
        end
end


function [xs, us, xs_bar, us_bar, t, modes, costs] = SimpleIntegrationMPC(obj, t0, tf, x0, x_star, usn, mpc_solver) %TODO
    local_mpc_solver = mpc_solver; % To avoid corrupting the mpc starting data
    num_integration_steps = ceil((tf - t0) / obj.h_step);
    assert(length(obj.real_states_map(1).x) == length(x0), 'x_state in Euler integration and real states has different size');
    number_of_variables = length(x0);
    xs = zeros(number_of_variables, number_of_integration_steps);
    xs(:, 1) = x0; 
    xs_bar = zeros(number_of_variables, number_of_integration_steps);
    t = zeros(number_of_integration_steps, 1);
    modes = zeros(1, number_of_integration_steps);
    costs = zeros(1, number_of_integration_steps);
    number_of_controllers = length(usn(t(1)));
    us = zeros(number_of_controllers, number_of_integration_steps);
    us_bar = zeros(number_of_controllers, number_of_integration_steps);
    for step = 1:number_of_integration_steps
        current_t = t(step);
        disp(current_t);
        current_x_state = xs(:, step);
        [us(:, step), modes(step), costs(step), local_mpc_solver] = local_mpc_solver.SolveMPC(current_t, x_star, usn, xs(:, step)); % Solve Controller MPC problem and get control input
        us_bar(:, step) = us(:, step) - usn(current_t);
        xs_bar(:, step) = xs(:, step) - x_star(current_t);
        if step < num_integration_steps
            delta_xs = obj.simulator.pointSimulator(xs(:, step), us(:, step));
            t(step + 1) = current_t + obj.h_step;
            xs(:, step + 1) = current_x_state + obj.h_step * delta_xs;  % + [normrnd(0,.005);normrnd(0,.005);normrnd(0,.005);0;0]; TODO: Implement noise if needed
        end
    end
end

% Euler-Integration for MPC
function [xs, xc, us, uc, xs_bar, us_bar, t, costs, fh_x_states, fh_considered_x] = IntegrateMPC(obj, t0, tf, xs0, xns, uns, mpc_solver) %TODO
    local_mpc_solver = mpc_solver; % To avoid corrupting the mpc starting data
    num_integration_steps = ceil((tf - t0) / obj.h_step);
    num_variables = length(xs0);
    xs = zeros(num_variables, num_integration_steps);
    xs(:, 1) = xs0; 
    xc = [];
    xs_bar = zeros(num_variables, num_integration_steps);
    t = zeros(num_integration_steps, 1);
    costs = zeros(1, num_integration_steps);
    number_of_controllers = length(ppval(uns, t(1)));
    us = zeros(number_of_controllers, num_integration_steps);
    uc = [];
    us_bar = zeros(number_of_controllers, num_integration_steps);
    fh_x_states = cell(num_integration_steps, 1);
    fh_considered_x = cell(num_integration_steps, 1);
    for step = 1:num_integration_steps
        disp(t(step));
        current_x_state = xs(:, step);
        [uc_i, us(:, step), costs(step), local_mpc_solver, fh_x_states{step}, fh_considered_x{step}] = local_mpc_solver.SolveMPC(t(step), xns, uns, xs(:, step)); % Solve Controller MPC problem and get control input
        fh_x = fh_x_states{step};
        xc = [xc, fh_x(:,1)];
        uc = [uc, uc_i];
        us_bar(:, step) = us(:, step) - ppval(uns, t(step));
        xs_bar(:, step) = xs(:, step) - ppval(xns, t(step));
        if step < num_integration_steps
            delta_xs = obj.simulator.pointSimulator(xs(:, step), us(:, step));
            t(step + 1) = t(step) + obj.h_step;
            xs(:, step + 1) = current_x_state + obj.h_step * delta_xs; 
        end
    end
end
end
    
end

