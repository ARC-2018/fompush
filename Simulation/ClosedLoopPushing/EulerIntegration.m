classdef EulerIntegration
    % Class to implement Euler Integration for certain desired
    % trajectories. It currently implements simple integration and MPC
    % solving
properties (Access = private)
    real_states_map; % Array with the real hybrid states to use to determine the dynamics of the Euler Integration step
    h_step;
    sim_count;
end
methods
function obj = EulerIntegration(real_states_map, h_step)
    assert(~isempty(real_states_map) > 0, 'Empty real_states array');
    obj.real_states_map = real_states_map;
    obj.h_step = h_step;
end

% function [] = GetScoresOneStep(obj, t0, x_star, u_star, mpc_solver)
%     local_mpc_solver = mpc_solver; % To avoid corrupting the mpc starting data
%     assert(length(obj.real_states_map(1).x) == length(x0), 'x_state in Euler integration and real states has different size');
%     x_state = zeros(length(x0), 1);
%     x_state(:, 1) = x0; 
%     number_of_controllers = length(u_star(t0));
%     u_bar = zeros(number_of_controllers, 1);
%         [modes, costs] = local_mpc_solver.SolveMPC(t0, x_star, u_star, x_state(:, step)); % Solve Controller MPC problem and get control input
%         for real_state_index = 1:(length(obj.real_states_map) + 1) % Determine which state this controller actually falls into
%             assert(real_state_index < length(obj.real_states_map) + 1, sprintf('At step %d of EulerIntegration, the real state to apply Euler integration cannot be determined', step));
%         end
% end


function [x_state, u_c_state, x_bar, u_bar, t, modes, costs] = SimpleIntegrationMPC(obj, t0, tf, x0, x_s_star, u_s_star, mpc_solver)
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
    number_of_controllers = length(u_s_star(t(1)));
    u_c_state = zeros(number_of_controllers, number_of_integration_steps);
    u_bar = zeros(number_of_controllers, number_of_integration_steps);
    for step = 1:number_of_integration_steps
        current_t = t(step);
        disp(current_t);
        current_x_state = x_state(:, step);
        [u_c_state(:, step), modes(step), costs(step), local_mpc_solver] = local_mpc_solver.SolveMPC(current_t, x_s_star, u_s_star, x_state(:, step)); % Solve Controller MPC problem and get control input
        u_bar(:, step) = u_c_state(:, step) - u_s_star(current_t);
        x_bar(:, step) = x_state(:, step) - x_s_star(current_t);
        for real_state_index = 1:(length(obj.real_states_map) + 1) % Determine which state this controller actually falls into
            assert(real_state_index < length(obj.real_states_map) + 1, sprintf('At step %d of EulerIntegration, the real state to apply Euler integration cannot be determined', step));
            if obj.real_states_map(real_state_index).CheckConstraints(x_state(:, step), u_c_state(:, step))
                disp(['Real state: ', obj.real_states_map(real_state_index).name]);
                delta_x = obj.real_states_map(real_state_index).GetMotionFunction(x_state(:, step), u_c_state(:, step));
                break;
            end
        end
        if step < number_of_integration_steps
            t(step + 1) = current_t + obj.h_step;
            x_state(:, step + 1) = current_x_state + obj.h_step * delta_x;  % + [normrnd(0,.005);normrnd(0,.005);normrnd(0,.005);0;0]; TODO: Implement noise if needed
        end
        %Save data 
        plot_data.xs_star_state(step,:) = x_s_star(t(step))';
        plot_data.us_star_state(step,:) = u_s_star(t(step));
        [xc_star, uc_star] = mpc_solver.simulator2controller(x_s_star(t(step)), u_s_star(t(step)));
        plot_data.xc_star_state(step,:) = xc_star';
        plot_data.uc_star_state(step,:) = uc_star';
        plot_data.us_state(step,:) = us;
        plot_data.uc_state(step,:) = u_c_state(:, step)';
    end
    %Post-Processing Data (Animation and plots)
    plot_data.t = t;
    plot_data.xs{2} = plot_data.xs_state;  
    plot_data.us{2} = plot_data.us_state; 
    plot_data.uc{2} = plot_data.uc_state; 
end

% Euler-Integration for MPC
function [x_s_state, u_s_state, x_bar, u_bar, t, modes, costs, animator] = IntegrateMPC(obj, t0, tf, x0_s, x_star, u_star, mpc_solver, animator)
    local_mpc_solver = mpc_solver; % To avoid corrupting the mpc starting data
    number_of_integration_steps = ceil((tf - t0) / obj.h_step);
    assert(obj.real_states_map.num_x == length(x0_s), 'x_state in Euler integration and real states has different size');
    num_x_s = length(x_star(0));
    x_s_state = zeros(num_x_s, number_of_integration_steps);
    x_s_state(:, 1) = x0_s; 
    x_bar = zeros(num_x_s, number_of_integration_steps);
    t = zeros(number_of_integration_steps, 1);
    modes = zeros(1, number_of_integration_steps);
    costs = zeros(1, number_of_integration_steps);
    num_u_s = length(u_star(0));
    u_s_state = zeros(num_u_s, number_of_integration_steps);
    u_bar = zeros(num_u_s, number_of_integration_steps);
    for step = 1:number_of_integration_steps
        current_t = t(step);
        disp(current_t);
        x_s = x_s_state(:, step);
        [x_c, ~] = mpc_solver.simulator2controller(x_s);
        [u_c, modes(step), costs(step), ~] = local_mpc_solver.SolveMPC(current_t, x_star, u_star, x_s_state(:, step)); % Solve Controller MPC problem and get control input
        [~, u_s_state(:, step)] = local_mpc_solver.controller2simulator(x_c, u_c);
        u_bar(:, step) = u_s_state(:, step) - u_star(current_t);
        x_bar(:, step) = x_s_state(:, step) - x_star(current_t);
        for real_state_index = 1:(length(obj.real_states_map) + 1) % Determine which state this controller actually falls into
            assert(real_state_index < length(obj.real_states_map) + 1, sprintf('At step %d of EulerIntegration, the real state to apply Euler integration cannot be determined', step));
%             if obj.real_states_map(real_state_index).CheckConstraints(x_state(:, step), u_state(:, step))
            if 1 == 1 % TODO: Add check if I can
                disp(['Real state: ', obj.real_states_map(real_state_index).name]);
                delta_x = obj.real_states_map(real_state_index).GetMotionFunction(x_s_state(:, step), u_s_state(:, step));
                break;
            end
        end
        if step < number_of_integration_steps
            t(step + 1) = current_t + obj.h_step;
            x_s_state(:, step + 1) = x_s + obj.h_step * delta_x;  % + [normrnd(0,.005);normrnd(0,.005);normrnd(0,.005);0;0]; TODO: Implement noise if needed
        end
        %Save data 
        animator.xs_star_state(step,:) = x_star(t(step))';
        animator.us_star_state(step,:) = u_star(t(step));
        [xc_star, uc_star] = mpc_solver.simulator2controller(x_star(t(step)), u_star(t(step))); %TODO: Review
        animator.xc_star_state(step,:) = xc_star';
        animator.uc_star_state(step,:) = uc_star';
        animator.us_state(step,:) = u_s_state(:, step)';
        animator.uc_state(step,:) = u_c';
    end
    %Post-Processing Data (Animation and plots)
    animator.t = t;
    animator.xs{2} = x_s_state';  
    animator.us{2} = u_s_state'; 
    animator.uc{2} = animator.uc_state; 
    animator.xs{1} = (x_s_state - x_bar)';  
    animator.us{1} = (u_s_state - u_bar)'; 
end
end
    
end

