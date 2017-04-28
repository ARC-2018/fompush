classdef EulerIntegration
    % Class to implement Euler Integration for certain desired
    % trajectories. It currently implements simple integration and MPC
    % solving
properties (Access = private)
    real_states_map; % Array with the real hybrid states to use to determine the dynamics of the Euler Integration step
    h_step;
end
methods
function obj = EulerIntegration(real_states_map, h_step)
    assert(~isempty(real_states_map) > 0, 'Empty real_states array');
    obj.real_states_map = real_states_map;
    obj.h_step = h_step;
end

function [] = GetScoresOneStep(obj, t0, x_star, u_star, mpc_solver)
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


function [x_state, u_state, x_bar, u_bar, t, modes, costs] = SimpleIntegrationMPC(obj, t0, tf, x0, x_star, u_star, mpc_solver)
    local_mpc_solver = mpc_solver; % To avoid corrupting the mpc starting data
    number_of_integration_steps = ceil((tf - t0) / obj.h_step);
    assert(length(obj.real_states_map(1).x) == length(x0), 'x_state in Euler integration and real states has different size');
    number_of_variables = length(x0);
    x_state = zeros(number_of_variables, number_of_integration_steps);
    x_state(:, 1) = x0; 
    x_bar = zeros(number_of_variables, number_of_integration_steps);
    t = zeros(number_of_integration_steps, 1);
    modes = zeros(1, number_of_integration_steps);
    costs = zeros(1, number_of_integration_steps);
    number_of_controllers = length(u_star(t(1)));
    u_state = zeros(number_of_controllers, number_of_integration_steps);
    u_bar = zeros(number_of_controllers, number_of_integration_steps);
    for step = 1:number_of_integration_steps
        current_t = t(step);
        disp(current_t);
        current_x_state = x_state(:, step);
        [u_state(:, step), modes(step), costs(step), local_mpc_solver] = local_mpc_solver.SolveMPC(current_t, x_star, u_star, x_state(:, step)); % Solve Controller MPC problem and get control input
        u_bar(:, step) = u_state(:, step) - u_star(current_t);
        x_bar(:, step) = x_state(:, step) - x_star(current_t);
        for real_state_index = 1:(length(obj.real_states_map) + 1) % Determine which state this controller actually falls into
            assert(real_state_index < length(obj.real_states_map) + 1, sprintf('At step %d of EulerIntegration, the real state to apply Euler integration cannot be determined', step));
            if obj.real_states_map(real_state_index).CheckConstraints(x_state(:, step), u_state(:, step))
                disp(['Real state: ', obj.real_states_map(real_state_index).name]);
                delta_x = obj.real_states_map(real_state_index).GetMotionFunction(x_state(:, step), u_state(:, step));
                break;
            end
        end
        if step < number_of_integration_steps
            t(step + 1) = current_t + obj.h_step;
            x_state(:, step + 1) = current_x_state + obj.h_step * delta_x;  % + [normrnd(0,.005);normrnd(0,.005);normrnd(0,.005);0;0]; TODO: Implement noise if needed
        end
    end
end

% Euler-Integration for MPC
function [x_state, u_state, x_bar, u_bar, t, modes, costs, fh_x_states, fh_considered_x] = IntegrateMPC(obj, t0, tf, x0, x_star, u_star, mpc_solver)
    local_mpc_solver = mpc_solver; % To avoid corrupting the mpc starting data
    number_of_integration_steps = ceil((tf - t0) / obj.h_step);
    assert(length(obj.real_states_map(1).x) == length(x0), 'x_state in Euler integration and real states has different size');
    number_of_variables = length(x0);
    x_state = zeros(number_of_variables, number_of_integration_steps);
    x_state(:, 1) = x0; 
    x_bar = zeros(number_of_variables, number_of_integration_steps);
    t = zeros(number_of_integration_steps, 1);
    modes = zeros(1, number_of_integration_steps);
    costs = zeros(1, number_of_integration_steps);
    number_of_controllers = length(u_star(t(1)));
    u_state = zeros(number_of_controllers, number_of_integration_steps);
    u_bar = zeros(number_of_controllers, number_of_integration_steps);
    u_bars = cell(number_of_integration_steps, 1);
    fh_x_states = cell(number_of_integration_steps, 1);
    fh_considered_x = cell(number_of_integration_steps, 1);
    for step = 1:number_of_integration_steps
        disp(t(step));
        current_x_state = x_state(:, step);
        [u_state(:, step), modes(step), costs(step), local_mpc_solver, fh_x_states{step}, fh_considered_x{step}] = local_mpc_solver.SolveMPC(t(step), x_star, u_star, x_state(:, step)); % Solve Controller MPC problem and get control input
%         [u_state(:, step), modes(step), costs(step), local_mpc_solver] = local_mpc_solver.SolveMPC(current_t, x_star, u_star, x_state(:, step)); % Solve Controller MPC problem and get control input
        u_bar(:, step) = u_state(:, step) - u_star(t(step));
        x_bar(:, step) = x_state(:, step) - x_star(t(step));
        if step < number_of_integration_steps
            for real_state_index = 1:(length(obj.real_states_map) + 1) % Determine which state this controller actually falls into
                assert(real_state_index < length(obj.real_states_map) + 1, sprintf('At step %d of EulerIntegration, the real state to apply Euler integration cannot be determined', step));
                if obj.real_states_map(real_state_index).CheckConstraints(x_state(:, step), u_state(:, step))
                    disp(['Real state: ', obj.real_states_map(real_state_index).name]);
                    delta_x = obj.real_states_map(real_state_index).GetMotionFunction(x_state(:, step), u_state(:, step));
                    break;
                end
            end
            t(step + 1) = t(step) + obj.h_step;
            x_state(:, step + 1) = current_x_state + obj.h_step * delta_x; % + [normrnd(0,.0005);normrnd(0,.0005);0;0]; %TODO: Implement noise if needed
        end
    end
end
end
    
end

