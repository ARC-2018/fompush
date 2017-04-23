classdef FOMSolver < MPCSolvers.MPCSolver
%FOMSOLVER Implementation of the MPCSolver that solves it by using the
%Family of Modes (FOM) approach.

properties (Access = private)
    hybrid_modes; % Array with the hybrid_modes to consider for the FOM.
    chameleon_mode;
    has_chameleon;
end

methods
function obj = FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, simulator2controller, controller2simulator, hybrid_modes, has_chameleon)
    obj = obj@MPCSolvers.MPCSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, simulator2controller, controller2simulator);
    assert(~isempty(hybrid_modes) > 0, 'Empty hybrid_modes array');
    obj.hybrid_modes = hybrid_modes;
    obj.has_chameleon = has_chameleon;
    obj.chameleon_mode = [];
end

%Function to only execute first step. Used for learning
function [f_values] = GetDataFirstStep(obj, current_t, u_star, current_x)
    number_of_modes = size(obj.hybrid_modes, 1);
    out_u_bar = cell(1, number_of_modes);
    f_values = Inf * ones(1, number_of_modes);
    for hybrid_mode_index = 1:number_of_modes
        hybrid_mode = current_hybrid_modes(hybrid_mode_index, :);
        optimization_problem = obj.GetOptimizationProblem(current_t, x_star, u_star, current_x, hybrid_mode);
        try
            [solved_optimization_problem, solvertime, f_values(hybrid_mode_index)] = optimization_problem.solve;
            out_u_bar{hybrid_mode_index} = solved_optimization_problem.vars.u.value;
        catch
            f_values(hybrid_mode_index) = Inf; % Disregard this solution
            fprintf('Opt. number %d not feasible\n', hybrid_mode_index);
        end
    end
    [min_cost, mode_index] = min(f_values);
    try
        u_bar = out_u_bar{mode_index}(1:obj.number_of_controllers, 1);
    catch
        error('Could not find a solution for any optimization problem');
    end
end

function [u_c, mode_index, min_cost, obj] = SolveMPC(obj, current_t, x_s_star, u_s_star, current_x_s)
    current_hybrid_modes = obj.hybrid_modes;
    if obj.has_chameleon
        current_hybrid_modes = [current_hybrid_modes; obj.chameleon_mode];
    end
    number_of_modes = size(current_hybrid_modes, 1);
    out_u_bar = cell(1, number_of_modes);
    out_x_bar = cell(1, number_of_modes);
    f_values = Inf * ones(1, number_of_modes);
    for hybrid_mode_index = 1:number_of_modes
        hybrid_mode = current_hybrid_modes(hybrid_mode_index, :);
        optimization_problem = obj.GetOptimizationProblem(current_t, x_s_star, u_s_star, current_x_s, hybrid_mode);
        try
            options = optimoptions('quadprog','Display','none');
            [solved_optimization_problem, solvertime, f_values(hybrid_mode_index)] = optimization_problem.solve;
            out_u_bar{hybrid_mode_index} = solved_optimization_problem.vars.u.value;
            out_x_bar{hybrid_mode_index} = solved_optimization_problem.vars.x.value;
        catch
            f_values(hybrid_mode_index) = Inf; % Disregard this solution
            fprintf('Opt. number %d not feasible\n', hybrid_mode_index);
        end
    end
    [min_cost, mode_index] = min(f_values);
    try
        u_c_bar = out_u_bar{mode_index}(:, 1);
%         x_bar = out_x_bar{hybrid_mode_index}(1:obj.number_of_variables, 1);
    catch
        error('Could not find a solution for any optimization problem');
    end
    disp([sprintf('Mode %d. First state: ', mode_index), obj.hybrid_states_map(current_hybrid_modes(mode_index, 1)).name]);
    [~, u_c_star] = obj.simulator2controller(x_s_star(current_t), u_s_star(current_t));
    u_c = u_c_bar + u_c_star;
    if obj.has_chameleon
        obj.chameleon_mode = [current_hybrid_modes(mode_index, 2:end), current_hybrid_modes(mode_index, 1)];
    end
end
end
methods
% methods (Access = private)
function [optimization_problem] = GetOptimizationProblem(obj, t0, x_s_star, u_star, x_s_0, hybrid_mode)
    number_of_steps = length(hybrid_mode); % Number of hybrid states in the hybrid_mode, it corresponds to the number of steps of the MPC problem
    [x_c_0, ~] = obj.simulator2controller(x_s_0);
    t = t0:obj.h_opt:(t0 + obj.h_opt * (number_of_steps - 1));
    u_lb = zeros(length(obj.u_lower_bound), number_of_steps);
    u_ub = zeros(length(obj.u_upper_bound), number_of_steps);
    for i = 1: number_of_steps
        [~, u_c] = obj.simulator2controller(x_s_star(t(i)), u_star(t(i)));
        u_lb(:,i) = -u_c - obj.u_lower_bound;
        u_ub(:,i) = -u_c + obj.u_upper_bound;
    end
    x_lb = - obj.x_lower_bound * ones(1, number_of_steps);
    x_ub = obj.x_upper_bound * ones(1, number_of_steps);
    optimization_problem = MixedIntegerConvexProgram(false); % Define optimization program
    % The arguments for the function are (name, type_, size_, lower_bound, upper_bound, start_value)
    optimization_problem = optimization_problem.addVariable('x', 'C', [obj.number_of_variables, number_of_steps], x_lb, x_ub);
    optimization_problem = optimization_problem.addVariable('u', 'C', [obj.number_of_controllers, number_of_steps], u_lb, u_ub);
    % Loop through steps of opt. program
    
    H = zeros(optimization_problem.nv, optimization_problem.nv);
    for step_index = 1:number_of_steps
        [x_c, u_c] = obj.simulator2controller(x_s_star(t(step_index)), u_star(t(step_index)));
        hybrid_state_index = hybrid_mode(step_index);
        hybrid_state = obj.hybrid_states_map(hybrid_state_index); % Todo change if it slows everything down
        [optimization_problem, H, Ai, bi, Ae, be] = hybrid_state.buildConstraints(optimization_problem, step_index, x_c_0, x_c, u_c, number_of_steps, obj.h_opt, obj.Q_MPC, obj.R_MPC, H);
        optimization_problem = optimization_problem.addLinearConstraints(Ai, bi, Ae, be);
    end
    optimization_problem = optimization_problem.addCost(H, [], []);
end

end
    
end

