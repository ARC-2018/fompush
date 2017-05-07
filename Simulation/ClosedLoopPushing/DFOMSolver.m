classdef DFOMSolver < MPCSolvers.MPCSolver
%FOMSOLVER Implementation of the MPCSolver that solves it by using the
%Family of Modes (FOM) approach.

properties (Access = private)
    hybrid_modes; % Array with the hybrid_modes to consider for the FOM.
    curvatures;
    best_modes;
    chameleon_mode;
    epsilon = 0.001;
end

methods
function obj = DFOMSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes, best_modes, curvatures)
    obj = obj@MPCSolvers.MPCSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt);
    assert(~isempty(hybrid_modes) > 0, 'Empty hybrid_modes array');
    obj.hybrid_modes = hybrid_modes;
    obj.chameleon_mode = [];
    obj.curvatures = curvatures;
    obj.best_modes = best_modes;
end

function [u_state, mode_index, min_cost, obj, chosen_x, considered_x] = SolveMPC(obj, current_t, x_star, u_star, current_x_state)
    options = optimoptions('quadprog','Display','none'); %TODO: Never used
    x_s = x_star([current_t - obj.epsilon, current_t, current_t + obj.epsilon]);
    x = x_s(1, :);
    y = x_s(2, :);
    dx = diff(x);
    dy = diff(y);
    ddx = diff(dx);
    ddy = diff(dy);
    num   = dx .* ddy - ddx .* dy;
    denom = dx .* dx + dy .* dy;
    denom = sqrt(denom);
    denom = denom .* denom .* denom;
    curvatur = num ./ denom;
    if (abs(curvatur) > 20)
        current_best_modes = obj.best_modes(:,end);
    else
        tmp = abs(obj.curvatures-curvatur(1));
        [~, idx] = min(tmp);
        current_best_modes = obj.best_modes(:,idx);
    end
    current_hybrid_modes = [obj.hybrid_modes(current_best_modes, :); obj.chameleon_mode];
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
    obj.chameleon_mode = [current_hybrid_modes(mode_index, 2:end), 1];
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
    t = [t, t0 + obj.h_opt * number_of_steps];
%     x_lb = -[100;100;100;100] * ones(1, number_of_steps);
%     x_ub = [100;100;100;100] * ones(1, number_of_steps);
    optimization_problem = MixedIntegerConvexProgram(false); % Define optimization program
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
            number_of_motion_cone_constraints = size(D,1);
            A_constraint = zeros(number_of_motion_cone_constraints, optimization_problem.nv);
        else
            [A, B, F, D, E, g] = hybrid_state.GetLinearMatrices(x_star(t(step_index)), u_star(t(step_index)));
            %% Dynamic Constraints
            A_bar = eye(size(A)) + obj.h_opt * A;
            B_bar = obj.h_opt * B;
            assert(size(A_bar, 1) == size(B_bar, 1), 'A_bar row number: %d and B_bar row number: %d mismatch', size(A_bar, 1), size(B_bar, 1));
            assert(size(A_bar, 1) == size(A_bar, 2), 'A_bar row number: %d and A_bar column number: %d mismatch', size(A_bar, 1), size(A_bar, 2));
            assert(size(A_bar, 1) == obj.number_of_variables, 'A_bar row number: %d and number of variables: %d mismatch', size(A_bar, 1), obj.number_of_variables);
            assert(size(B_bar, 2) == obj.number_of_controllers, 'B_bar column number: %d and number of controllers: %d mismatch', size(B_bar, 2), obj.number_of_controllers);
            assert(size(E, 1) == size(D, 1), 'E row number: %d, and D row number: %d mismatch', size(E, 1), size(D, 1));
            assert(size(D, 1) == size(g, 1), 'D row number: %d, and g row number: %d mismatch', size(D, 1), size(g, 1));
            assert(size(E, 2) == size(A_bar, 2), 'E column number: %d, and A_bar column number: %d mismatch', size(E, 2), size(A_bar, 2));
            assert(size(D, 2) == size(B_bar, 2), 'D column number: %d, and B_bar column number: %d mismatch', size(D, 2), size(B_bar, 2));
            A_motion(:,optimization_problem.vars.x.i(1:obj.number_of_variables, step_index-1)) = -A_bar;
            A_motion(:,optimization_problem.vars.u.i(1:obj.number_of_controllers, step_index)) = -B_bar;
%             b_motion = zeros(size(A_motion, 1),1);
            b_motion = x_star(t(step_index)) - x_star(t(step_index + 1)) + obj.h_opt * F;
            number_of_motion_cone_constraints = size(D,1);
            A_constraint = zeros(number_of_motion_cone_constraints, optimization_problem.nv);
            A_constraint(:, optimization_problem.vars.x.i(1:obj.number_of_variables, step_index)) = E;
        end
        %Final Cost
        if step_index == number_of_steps
%             H(optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index), optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index)) = 10 * dare(A, B, obj.Q_MPC, obj.R_MPC);
            H(optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index), optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index)) = obj.Q_MPC_final;
        end
        optimization_problem = optimization_problem.addCost(H, [], []);
        A_constraint(:, optimization_problem.vars.u.i(1:obj.number_of_controllers, step_index)) = D;
        b_constraint = g;
        optimization_problem = optimization_problem.addLinearConstraints(A_constraint, b_constraint, A_motion, b_motion);
    end
end

end
    
end

