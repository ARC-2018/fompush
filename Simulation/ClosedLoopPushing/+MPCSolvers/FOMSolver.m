classdef FOMSolver < MPCSolvers.MPCSolver
%FOMSOLVER Implementation of the MPCSolver that solves it by using the
%Family of Modes (FOM) approach.

properties (Access = private)
    hybrid_modes; % Array with the hybrid_modes to consider for the FOM.
    chameleon_mode;
    has_chameleon;
end

methods
function obj = FOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes, has_chameleon)
    obj = obj@MPCSolvers.MPCSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt);
    assert(~isempty(hybrid_modes) > 0, 'Empty hybrid_modes array');
    obj.hybrid_modes = hybrid_modes;
    obj.has_chameleon = has_chameleon;
    obj.chameleon_mode = [];
end

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

function [u_state, mode_index, min_cost, obj, chosen_x, considered_x] = SolveMPC(obj, current_t, x_star, u_star, current_x_state)
    options = optimoptions('quadprog','Display','none'); %TODO: Never used
    current_hybrid_modes = obj.hybrid_modes;
    if obj.has_chameleon
        current_hybrid_modes = [current_hybrid_modes; obj.chameleon_mode];
    end
    number_of_modes = size(current_hybrid_modes, 1);
    out_u_bar = cell(1, number_of_modes);
    considered_x = cell(1, number_of_modes);
    f_values = Inf * ones(1, number_of_modes);
    t = current_t:obj.h_opt:(current_t + obj.h_opt * length(current_hybrid_modes(1, :)));
    for hybrid_mode_index = 1:number_of_modes
        hybrid_mode = current_hybrid_modes(hybrid_mode_index, :);
        optimization_problem = obj.GetOptimizationProblem(current_t, x_star, u_star, current_x_state, hybrid_mode);
        try
            [solved_optimization_problem, solvertime, f_values(hybrid_mode_index)] = optimization_problem.solve;
            out_u_bar{hybrid_mode_index} = solved_optimization_problem.vars.u.value;
            considered_x{hybrid_mode_index} = [current_x_state, solved_optimization_problem.vars.x.value + x_star(t(2:end))];
        catch
            f_values(hybrid_mode_index) = Inf; % Disregard this solution
            fprintf('Opt. number %d not feasible\n', hybrid_mode_index);
        end
    end
    [min_cost, mode_index] = min(f_values);
%     x_values = problems{mode_index}.vars.x.value
%     u_values = problems{mode_index}.vars.u.value
%     f_values
    try
        u_bar = out_u_bar{mode_index}(1:obj.number_of_controllers, 1);
    catch
        error('Could not find a solution for any optimization problem');
    end
    disp([sprintf('Mode %d. First state: ', mode_index), obj.hybrid_states_map(current_hybrid_modes(mode_index, 1)).name]);
    u_state = u_bar + u_star(current_t);
    chosen_x = considered_x{mode_index};
    if obj.has_chameleon
        obj.chameleon_mode = [current_hybrid_modes(mode_index, 2:end), 1];
    end
end
end
methods
% methods (Access = private)
function [optimization_problem] = GetOptimizationProblem(obj, t0, x_star, u_star, x_0_state, hybrid_mode)
    number_of_steps = length(hybrid_mode); % Number of hybrid states in the hybrid_mode, it corresponds to the number of steps of the MPC problem
    t = t0:obj.h_opt:(t0 + obj.h_opt * (number_of_steps - 1));
    u_lb = -u_star(t) - obj.u_lower_bound * ones(1, number_of_steps);
    u_ub = -u_star(t) + obj.u_upper_bound * ones(1, number_of_steps);
    x_lb = [0;0;0;1] * ones(1, number_of_steps) .* -x_star(t) - obj.x_lower_bound * ones(1, number_of_steps);
    x_ub = [0;0;0;1] * ones(1, number_of_steps) .* -x_star(t) + obj.x_upper_bound * ones(1, number_of_steps);
    x_s = x_star(t);
%     x_lb = -[100;100;100;100] * ones(1, number_of_steps);
%     x_ub = [100;100;100;100] * ones(1, number_of_steps);
    optimization_problem = MixedIntegerConvexProgram(false);
    t = [t, t0 + obj.h_opt * number_of_steps]; % Define optimization program
    % The arguments for the function are (name, type_, size_, lower_bound, upper_bound, start_value)
    optimization_problem = optimization_problem.addVariable('x', 'C', [obj.number_of_variables, number_of_steps], x_lb, x_ub);
    optimization_problem = optimization_problem.addVariable('u', 'C', [obj.number_of_controllers, number_of_steps], u_lb, u_ub);
    % Loop through steps of opt. program
    for step_index = 1:number_of_steps;
        hybrid_state_index = hybrid_mode(step_index);
        hybrid_state = obj.hybrid_states_map(hybrid_state_index); % Todo change if it slows everything down
        %% Cost
        H = zeros(optimization_problem.nv, optimization_problem.nv);
        H(optimization_problem.vars.x.i(1:length(obj.Q_MPC), step_index), optimization_problem.vars.x.i(1:length(obj.Q_MPC), step_index)) = obj.Q_MPC;
        H(optimization_problem.vars.u.i(1:length(obj.R_MPC), step_index), optimization_problem.vars.u.i(1:length(obj.R_MPC), step_index)) = obj.R_MPC;
        A_motion = zeros(obj.number_of_variables, optimization_problem.nv);
        A_motion(:,optimization_problem.vars.x.i(1:obj.number_of_variables, step_index)) = eye(obj.number_of_variables);
        % TODO: Should add constraint to bound the first value to it's
        % real value?
        if step_index == 1
            delta_x0 = x_0_state - x_star(t(1));
            [B, D, g] = hybrid_state.GetInitialStateMatrices(x_0_state, u_star(t(1)));
            %% Add nonlinear dynamic constraint
            A_motion(:,optimization_problem.vars.u.i(1:obj.number_of_controllers, 1)) = -obj.h_opt * B;
            b_motion = obj.h_opt * B * u_star(t(1)) + delta_x0 + x_star(t(1)) - x_star(t(2));
%             [~, ~, F_star, ~, ~, ~] = hybrid_state.GetLinearMatrices(x_star(t(1)), u_star(t(1)));
%             b_motion = obj.h_opt * B * u_star(t(1)) -obj.h_opt * F_star + delta_x0;
            number_of_motion_cone_constraints = size(D,1);
            A_constraint = zeros(number_of_motion_cone_constraints, optimization_problem.nv);
        else
            [A, B, F, D, E, g] = hybrid_state.GetLinearMatrices(x_star(t(step_index)), u_star(t(step_index)));
            %% Dynamic Constraints
            A_bar = eye(size(A)) + obj.h_opt * A;
            B_bar = obj.h_opt * B;
            A_motion(:,optimization_problem.vars.x.i(1:obj.number_of_variables, step_index-1)) = -A_bar;
            A_motion(:,optimization_problem.vars.u.i(1:obj.number_of_controllers, step_index)) = -B_bar;
            b_motion = x_star(t(step_index)) - x_star(t(step_index + 1)) + obj.h_opt * F;
%             b_motion = zeros(size(A_motion, 1),1);
            number_of_motion_cone_constraints = size(D,1);
            A_constraint = zeros(number_of_motion_cone_constraints, optimization_problem.nv);
            A_constraint(:, optimization_problem.vars.x.i(1:obj.number_of_variables, step_index)) = E;
        end
        %Final Cost
        if step_index == number_of_steps
            H(optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index), optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index)) = 10 * dare(A, B, obj.Q_MPC, obj.R_MPC);
%             H(optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index), optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index)) = obj.Q_MPC_final;
        end
        optimization_problem = optimization_problem.addCost(H, [], []);
        A_constraint(:, optimization_problem.vars.u.i(1:obj.number_of_controllers, step_index)) = D;
        b_constraint = g;
        optimization_problem = optimization_problem.addLinearConstraints(A_constraint, b_constraint, A_motion, b_motion);
    end
end

end
    
end

