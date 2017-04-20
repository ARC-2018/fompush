classdef MPCSolver
%MPCPROBLEM Interface for the MPC Problem solver

properties
    hybrid_states_map; % It contains the hybrid states to consider in the MPC
    number_of_variables;
    number_of_controllers;
    simulator2controller;
    controller2simulator;
end

properties(Access = protected)
    Q_MPC; % Cost matrix for the instantaneous x_state cost of the MPC.
    Q_MPC_final; % Cost matrix for the x_state cost of the final state of the MPC.
    R_MPC; % Cost matrix for the instantaneous u_state cost of the MPC.
    u_lower_bound;
    u_upper_bound;
    x_lower_bound;
    x_upper_bound;
    h_opt;
end

methods
function obj = MPCSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, simulator2controller, controller2simulator)
    obj.hybrid_states_map = hybrid_states_map;
    obj.number_of_variables = hybrid_states_map(1).num_x;
    obj.number_of_controllers = hybrid_states_map(1).num_u;
    obj.Q_MPC = Q_MPC;
    obj.Q_MPC_final = Q_MPC_final;
    obj.R_MPC = R_MPC;
    obj.u_lower_bound = u_lower_bound;
    obj.u_upper_bound = u_upper_bound;
    obj.x_lower_bound = x_lower_bound;
    obj.x_upper_bound = x_upper_bound;
    obj.h_opt = h_opt;
    obj.simulator2controller = simulator2controller;
    obj.controller2simulator = controller2simulator;
end
end
methods (Abstract)
    [u_state, mode_index, min_cost] = SolveMPC(obj, current_t, x_star, u_star, current_x);
end
    
end
