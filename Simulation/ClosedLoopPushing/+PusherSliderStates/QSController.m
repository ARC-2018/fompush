classdef QSController
%Groups all the controller properties

properties
    num_x = 4;
    num_y = 2;
end

methods
    function optimization_problem = buildConstraints(obj, optimization_problem, index, x0, x_c, u_c, number_of_steps, h_opt, Q_MPC, R_MPC)
        % TODO: Should add constraint to bound the first value to it's real value?
        A_motion = zeros(obj.num_x, optimization_problem.nv);
        A_motion(:,optimization_problem.vars.x.i(1:obj.num_x, step_index)) = eye(obj.num_x);
        if index == 1
            delta_x0 = x0 - x_c;
            [B, F, D, g] = obj.GetInitialStateMatrices(x0, x_c, u_c);
            %% Add nonlinear dynamic constraint
            F_bar = delta_x0 + h_opt * F;
            B_bar = h_opt * B;
            %Add constraint (Modify existing dummy constraint)
            A_motion(:,optimization_problem.vars.u.i(1:obj.num_u, 1)) = -B_bar;
            b_motion = F_bar;
            number_of_motion_cone_constraints = size(D,1);
            A_constraint = zeros(number_of_motion_cone_constraints, optimization_problem.nv);
        else
            [A, B, D, E, g] = obj.GetLinearMatrices(x_c, u_c);
            %% Dynamic Constraints
            A_bar = eye(size(A)) + h_opt * A;
            B_bar = h_opt * B;
            A_motion(:,optimization_problem.vars.x.i(1:obj.num_x, index-1)) = -A_bar;
            A_motion(:,optimization_problem.vars.u.i(1:obj.num_u, index)) = -B_bar;
            b_motion = zeros(size(A_motion, 1),1);
            number_of_motion_cone_constraints = size(D,1);
            A_constraint = zeros(number_of_motion_cone_constraints, optimization_problem.nv);
            A_constraint(:, optimization_problem.vars.x.i(1:obj.num_x, index)) = E;
        end 
        A_constraint(:, optimization_problem.vars.u.i(1:obj.num_u, index)) = D;
        b_constraint = g;
        optimization_problem = optimization_problem.addLinearConstraints(A_constraint, b_constraint, A_motion, b_motion);
        %% Cost
        H = zeros(optimization_problem.nv, optimization_problem.nv);
        H(optimization_problem.vars.x.i(1:length(Q_MPC), step_index), optimization_problem.vars.x.i(1:length(Q_MPC), step_index)) = Q_MPC;
        H(optimization_problem.vars.u.i(1:length(R_MPC), step_index), optimization_problem.vars.u.i(1:length(R_MPC), step_index)) = R_MPC;
        %Final Cost
        if step_index == number_of_steps
            Q_final = 10 * dare(A_bar, B_bar, Q_MPC, R_MPC);
            H(optimization_problem.vars.x.i(1:obj.num_x, step_index), optimization_problem.vars.x.i(1:obj.num_x, step_index)) = Q_final;
        end
        optimization_problem = optimization_problem.addCost(H, [], []);
    end
end
methods (Abstract)
    [B, F, D, g] = GetInitialStateMatrices(obj, x0, x_star, u_star);
    [A, B, D, E, g] = GetLinearMatrices(obj, x_star, u_star);
end
    
end

