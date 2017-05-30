classdef Force1pointSolver  < MPCSolvers.MPCSolver
%FOMSOLVER Implementation of the MPCSolver that solves it by using the
%Family of Modes (FOM) approach.

properties (Access = private)
    hybrid_modes; % Array with the hybrid_modes to consider for the FOM.
    chameleon_mode;
    has_chameleon;
    lp_functions;
end

methods
function obj = Force1pointSolver(Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes)
    obj = obj@MPCSolvers.MPCSolver([struct('x', [1 1 1 1], 'u', [1 1 1])], Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt);
    assert(~isempty(hybrid_modes) > 0, 'Empty hybrid_modes array');
    obj.hybrid_modes = hybrid_modes;
    obj.chameleon_mode = [];
    obj.lp_functions = Oneptpusher.Controller();
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

function [uc, us, min_cost, obj, chosen_xc, considered_xc] = SolveMPC(obj, current_t, xns, uns, current_xs)
    options = optimoptions('quadprog','Display','none'); %TODO: Never used
    current_hybrid_modes = obj.hybrid_modes;
    current_hybrid_modes = [current_hybrid_modes; obj.chameleon_mode];
    number_of_modes = size(current_hybrid_modes, 1);
    out_uc_bar = cell(1, number_of_modes);
    considered_xc = cell(1, number_of_modes);
    f_values = Inf * ones(1, number_of_modes);
    t = current_t:obj.h_opt:(current_t + obj.h_opt * length(current_hybrid_modes(1, :)));
    for hybrid_mode_index = 1:number_of_modes
        hybrid_mode = current_hybrid_modes(hybrid_mode_index, :);
        op = obj.GetOptimizationProblem(current_t, xns, uns, current_xs, hybrid_mode);
        try
            [solved_optimization_problem, solvertime, f_values(hybrid_mode_index)] = op.solve;
            out_uc_bar{hybrid_mode_index} = solved_optimization_problem.vars.u.value;
            considered_xc{hybrid_mode_index} = [obj.lp_functions.coordinateTransformSC(current_xs), solved_optimization_problem.vars.x.value + obj.GetXC(xns, t(2:end))];
        catch ME
            display(ME.identifier);
            f_values(hybrid_mode_index) = Inf; % Disregard this solution
            fprintf('Opt. number %d not feasible\n', hybrid_mode_index);
        end
    end
    [min_cost, mode_index] = min(f_values);
    uc_bar = out_uc_bar{mode_index}(1:obj.number_of_controllers, 1);
    disp(sprintf('Mode %d. First state: ', mode_index));
    chosen_xc = considered_xc{mode_index};
    [~, unc] = obj.lp_functions.fullTransformSC(ppval(xns, current_t), ppval(uns, current_t));
    uc = uc_bar + unc;
    us = obj.lp_functions.force2Velocity(chosen_xc(:, 1), uc);
    obj.chameleon_mode = [current_hybrid_modes(mode_index, 2:end), 1];
end

function xnc = GetXC(obj, xns, t)
    xnc = [];
    for ti = t
        xnc = [xnc, obj.lp_functions.coordinateTransformSC(ppval(xns, ti))];
    end
end
    
function [xnc, unc] = TransformTrajectories(obj, xns, uns, t)
    xnc=[];
    unc=[];
    for ti = t
        [xnci, unci] = obj.lp_functions.fullTransformSC(ppval(xns, ti), ppval(uns, ti));
        xnc = [xnc, xnci];
        unc = [unc, unci];
    end
end
    
function [op] = GetOptimizationProblem(obj, t0, xns, uns, x0s, hybrid_mode)
    num_steps = length(hybrid_mode); % Number of hybrid states in the hybrid_mode, it corresponds to the number of steps of the MPC problem
    t = t0:obj.h_opt:(t0 + obj.h_opt * num_steps);
    [xnc, unc] = obj.TransformTrajectories(xns, uns, t);
%     u_lb = -unc(:, 1:(end - 1)) - obj.u_lower_bound * ones(1, num_steps);
%     u_ub = -unc(:, 1:(end - 1)) + obj.u_upper_bound * ones(1, num_steps);
%     x_lb = [0;0;0;1] * ones(1, num_steps) .* -xnc(:, 2:end) - obj.x_lower_bound * ones(1, num_steps);
%     x_ub = [0;0;0;1] * ones(1, num_steps) .* -xnc(:, 2:end) + obj.x_upper_bound * ones(1, num_steps);
    u_lb = -1000 * ones(3, num_steps);
    u_ub = 1000 * ones(3, num_steps);
    x_lb = -1000 * ones(4, num_steps);
    x_ub = 1000 * ones(4, num_steps);
    op = MixedIntegerConvexProgram(false); % Define optimization program
    op = op.addVariable('x', 'C', [obj.number_of_variables, num_steps], x_lb, x_ub); % (name, type_, size_, lower_bound, upper_bound, start_value)
    op = op.addVariable('u', 'C', [obj.number_of_controllers, num_steps], u_lb, u_ub);
    
    for step_index = 1:num_steps; % Loop through steps of opt. program
        mode_index = hybrid_mode(step_index);
        H = zeros(op.nv, op.nv); % Cost
        H(op.vars.x.i(1:length(obj.Q_MPC), step_index), op.vars.x.i(1:length(obj.Q_MPC), step_index)) = obj.Q_MPC;
        H(op.vars.u.i(1:length(obj.R_MPC), step_index), op.vars.u.i(1:length(obj.R_MPC), step_index)) = obj.R_MPC;
        
        %% Motion
        A_motion = zeros(obj.number_of_variables, op.nv);
        A_motion(:,op.vars.x.i(1:obj.number_of_variables, step_index)) = eye(obj.number_of_variables);
        b_motion = zeros(obj.number_of_variables, 1);
        [A, B] = obj.lp_functions.motionConstraintMatrices(xnc(:, step_index), unc(:, step_index));
        A_bar = eye(size(A)) + obj.h_opt * A;
        B_bar = obj.h_opt * B;
        %Special case of initial conditions
        if step_index ~=1
            A_motion(:,op.vars.x.i(1:obj.number_of_variables, step_index-1))= -A_bar;
        else
            delta_xc = obj.lp_functions.coordinateTransformSC(x0s) - xnc(1);
            b_motion = A_bar * delta_xc;
        end
        A_motion(:, op.vars.x.i(1:obj.number_of_variables, step_index))  = eye(obj.number_of_variables);
        A_motion(:, op.vars.u.i(1:obj.number_of_controllers, step_index))=  -B_bar;
        
        %% Force
        [Aeq, beq, Ain, bin] = obj.lp_functions.forceConstraintMatrices(xnc(:, step_index), unc(:, step_index), mode_index);
        Ain_force = zeros(size(Ain,1), op.nv);
        Ain_force(:,op.vars.u.i(:,step_index)) = Ain;
        bin_force = bin;
        
        Aeq_force = zeros(size(Aeq,1), op.nv);
        Aeq_force(:,op.vars.u.i(:,step_index)) = Aeq;
        beq_force = beq;
        
        %Final Cost
        if step_index == num_steps
            H(op.vars.x.i(1:length(obj.Q_MPC_final), step_index), op.vars.x.i(1:length(obj.Q_MPC_final), step_index)) = 10 * dare(A, B, obj.Q_MPC, obj.R_MPC);
%             H(optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index), optimization_problem.vars.x.i(1:length(obj.Q_MPC_final), step_index)) = obj.Q_MPC_final;
        end
        op = op.addCost(H, [], []);
        op = op.addLinearConstraints(Ain_force, bin_force, [A_motion; Aeq_force], [b_motion; beq_force]);
    end
end

end
    
end