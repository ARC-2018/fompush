classdef LPSolver < MPCSolvers.MPCSolver
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
function obj = LPSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt, hybrid_modes, best_modes, curvatures)
    obj = obj@MPCSolvers.MPCSolver(hybrid_states_map, Q_MPC, Q_MPC_final, R_MPC, u_lower_bound, u_upper_bound, x_lower_bound, x_upper_bound, h_opt);
    assert(~isempty(hybrid_modes) > 0, 'Empty hybrid_modes array');
    obj.hybrid_modes = hybrid_modes;
    obj.chameleon_mode = [];
    obj.curvatures = curvatures;
    obj.best_modes = best_modes;
    
    
            obj.symbolicLinearize;
            obj.buildNominalTrajectory;
end


function obj = symbolicLinearize(obj)
    % Builds A and B linear matrices
    %
    %   Parameters:
    %   A       -   Continuous linear state matrix
    %   B       -   Continuous linear input matrix
    %   A_bar   -   Discrete linear state matrix
    %   B_bar   -   Discrete linear input matrix
    %   A,B:    -   dx_bar = A*x_bar + B*u_bar
    %   A_bar,B_bar:    -   x_bar(i+1) = A_bar*x_bar(i) + B_bar*u_bar(i)
    %% Symbolic variables
    syms x y theta a xp yp ry
    syms fx fy m fn1 fn2 ft1 ft2 ry_dot
    %% Build states       
    xc = [x;y;theta;ry];
    uc = [fn1;fn2;ft1; ft2; ry_dot];
    fn = [fn1;fn2];
    ft = [ft1; ft2];
    %% DCM Matrices
    Cbi = Helper.C3_2d(theta);
    %% Kinematics
    sign_vec = [1 -1]*1;
    rx = -obj.a/2;
    for lv1=1:2 %lv1 represents contact point 1 and 2
        rb{lv1} = [rx*1;ry*1]+sign_vec(lv1)*[0;obj.d*1]; %position of contact point lv1
        Jb{lv1} = [1 0 -rb{lv1}(2);... 
                     0 1 rb{lv1}(1)];
        for lv2=1:2 %lv2 represents left of right border (FC/MC)
            n{lv1} = [1;0];
            t{lv1} = [0;1];
            N{lv1} = transpose(n{lv1})*Jb{lv1};
            T{lv1} = transpose(t{lv1})*Jb{lv1};
        end
    end
    %Motion equations (nonlinear)
    N_tot = [N{1};N{2}];
    T_tot = [T{1};T{2}];
    Vb = obj.A_ls*(transpose(N_tot)*fn + transpose(T_tot)*ft );
    C_tilde = [Cbi 0;0 0 1];
    f_non1 = C_tilde*Vb;
    f_non2 = ry_dot;
    f_non = [f_non1;f_non2];
    %Linearization
    A = jacobian(f_non,xc);
    B = jacobian(f_non,uc);
    % Substitute equilibrium states
    A = subs(A,{x,y,theta,ry},{obj.xc_eq(1),obj.xc_eq(2),obj.xc_eq(3), obj.xc_eq(4)});
    B = subs(B,{x,y,theta,ry},{obj.xc_eq(1),obj.xc_eq(2),obj.xc_eq(3), obj.xc_eq(4)});
    A = subs(A,{fn1,fn2,ft1,ft2,ry_dot},{obj.uc_eq(1),obj.uc_eq(2),obj.uc_eq(3),obj.uc_eq(4),obj.uc_eq(5)});
    B = subs(B,{fn1,fn2,ft1,ft2,ry_dot},{obj.uc_eq(1),obj.uc_eq(2),obj.uc_eq(3),obj.uc_eq(4),obj.uc_eq(5)});
    %Convert to double type
    A=double(A);
    B=double(B);
    %Set properties
    obj.A_linear = double(A);
    obj.B_linear = double(B);
    obj.A_bar = double(eye(4)+obj.h*A);
    obj.B_bar = double(obj.h*B);
end

%solve MPC optimization problem
        function [uc] = solveMPC(obj, xc, t)
            % Returns controller input (force control) given initial
            % condition and time index
            %
            %   USAGE
            %   [uc] = solveMPC(obj, xc, t)
            %
            %   Parameters:
            %   t         -   1x1 simulation time index
            %   xc        -   4x1 controller state [x;y;theta;ry]. (pose of slider AND relative pose of center of pusher)
            
            %Define variables
            x = xc(1);
            y = xc(2);
            theta = xc(3);
            ry = xc(4);
            %Nominal coordinates
            [xcStar, ucStar, ~, ~] = obj.getStateNominal(t);
            %Build error state vector
            delta_xc = xc - xcStar;
            %Loop through family of modes
            fVal = [];
            for lv2=1:3
            %Build optimization program
            Opt{lv2} = obj.buildProgram();
                %Loop through steps of MPC
                for lv1=1:obj.steps
                    %Add cost
                    Opt{lv2} = obj.buildCost(Opt{lv2}, lv1);
                    %Add dynamic constraints
                    Opt{lv2} = obj.buildDynConstraints(Opt{lv2}, lv1);
                    %Add mode independent constraints
                    Opt{lv2} = obj.buildModeIndepConstraints(Opt{lv2}, lv1);
                    if lv1==1
                        if lv2==1 
                            %Add mode dependant constraints
                            Opt{lv2} = obj.buildModeDepConstraints(Opt{lv2}, lv1, 1);%1: Sticking, 2: Sliding left(up), 3: Slide right(down)
                        elseif lv2==2 
                            %Add mode dependant constraints
                            Opt{lv2} = obj.buildModeDepConstraints(Opt{lv2}, lv1, 2);%1: Sticking, 2: Sliding left(up), 3: Slide right(down)
                        elseif lv2==3 
                            %Add mode dependant constraints
                            Opt{lv2} = obj.buildModeDepConstraints(Opt{lv2}, lv1, 3);%1: Sticking, 2: Sliding left(up), 3: Slide right(down)
                        end
                    else
                        %Add mode dependant constraints
                        Opt{lv2} = obj.buildModeDepConstraints(Opt{lv2}, lv1, 1);%1: Sticking, 2: Sliding left(up), 3: Slide right(down)
                    end  
                end
                %Update initial conditions
                Opt{lv2}.beq(1:4) = zeros(4,1);
                Opt{lv2}.beq(1:4) = [obj.A_bar*delta_xc];
                % Solve Opt Program   
                options = optimoptions('quadprog','Display','none');
                [Opt{lv2}, solvertime{lv2}, fval{lv2}] = Opt{lv2}.solve;
                out_delta_u{lv2} = Opt{lv2}.vars.u.value';
                out_delta_x{lv2} = Opt{lv2}.vars.x.value';
                fVal = [fVal; fval{lv2}];
            end
            %Find mode schedule with lowest cost
            [minFOM indexFOM] = min(fVal);
            %Return first element of control sequence
            delta_u = out_delta_u{indexFOM}(1,1:5)';
            %Add feedforward and feedback controls together
            uc = delta_u + ucStar;
        end

function [u_state, mode_index, min_cost, obj, chosen_x, considered_x] = SolveMPC(obj, t_0, x_s_n, u_s_n, x_s_0)
    options = optimoptions('quadprog','Display','none'); %TODO: Never used
    current_hybrid_modes = obj.hybrid_modes;
    current_hybrid_modes = [current_hybrid_modes; obj.chameleon_mode];
    
    number_of_modes = size(current_hybrid_modes, 1);
    out_u_bar = cell(1, number_of_modes);
    considered_x = cell(1, number_of_modes);
    f_values = Inf * ones(1, number_of_modes);
    t = t_0:obj.h_opt:(t_0 + obj.h_opt * length(current_hybrid_modes(1, :)));
    for hybrid_mode_index = 1:number_of_modes
        hybrid_mode = current_hybrid_modes(hybrid_mode_index, :);
        optimization_problem = obj.GetOptimizationProblem(t_0, x_s_n, u_s_n, x_s_0, hybrid_mode);
        try
            [solved_optimization_problem, solvertime, f_values(hybrid_mode_index)] = optimization_problem.solve;
            out_u_bar{hybrid_mode_index} = solved_optimization_problem.vars.u.value;
            considered_x{hybrid_mode_index} = [x_s_0, solved_optimization_problem.vars.x.value + x_s_n(t(2:end))];
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
    disp([sprintf('Mode %d. First state: ', mode_index), obj.hybrid_states_map(current_hybrid_modes(mode_index, 1)).name]);
    u_state = u_bar + u_s_n(t_0);
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

function [o_nc, o_ui, o_bl, o_cl, o_br, o_cr] = InitNUBC(nc, ui, bl, cl, br, cr)
    o_nc = nc; o_ui = ui; o_bl = bl; o_cl = cl; o_br = br; o_cr = cr;
end

%Build mode dependant constraints
function Opt = buildModeDepConstraints(obj, Opt, lv1, mode)
    if mode==1
        [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(1,5,1,0,0,0); % Sticking Constraint
        Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '==', 'star', Opt);
        clear Bleft cLeft Bright cRight
    elseif mode==2
        [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(1,5,1,0,0,0); % Sliding left constraint
        Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '>', 'star', Opt);
        clear Bleft cLeft Bright cRight
        [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(1,[1,3],[0,1],0,[obj.nu_p, 0],0); % friction cone edge constraints
        Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '==', 'star', Opt);
        clear Bleft cLeft Bright cRight
        [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(1,[2,4],[0 1],0,[obj.nu_p 0],0); % friction cone edge constraints
        Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '==', 'star', Opt);
        clear Bleft cLeft Bright cRight
    elseif mode==3
        [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(1,5,1,0,0,0); % Sliding left constraint
        Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '<=', 'star', Opt);
        clear Bleft cLeft Bright cRight
        [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(1,[1,3],[0,1],0,[-obj.nu_p, 0],0); % friction cone edge constraints
        Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '==', 'star', Opt);
        clear Bleft cLeft Bright cRight 
        [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(1,[2,4],[0 1],0,[-obj.nu_p 0],0); % friction cone edge constraints
        Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '==', 'star', Opt);
        clear Bleft cLeft Bright cRight
    end
end
%Build dynamic constraints
function Opt = buildDynConstraints(obj, Opt, lv1)
    numConstraints = obj.num_xcStates;
    Aeq = zeros(numConstraints, Opt.nv);
    %Special case of initial conditions
    if lv1 ~=1
    Aeq(:,Opt.vars.x.i(1:obj.num_xcStates,lv1-1))= -obj.A_bar;
    end
    Aeq(:,Opt.vars.x.i(1:obj.num_xcStates,lv1))  = eye(obj.num_xcStates);
    Aeq(:,Opt.vars.u.i(1:obj.num_ucStates,lv1))=  -obj.B_bar;
    beq = zeros(obj.num_xcStates,1);
    Opt = Opt.addLinearConstraints([], [], Aeq, beq);
end
%Build mode independant constraints
function Opt = buildModeIndepConstraints(obj, Opt, lv1)

    [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(2,[1, 2],[],[],eye(2),[]); % positive normal force(s)
    Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '<=', 'star', Opt);
    clear Bleft cLeft Bright cRight
    [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(1,[1, 3],[0 1],[],[obj.nu_p 0],[]); % friction cone clamping of frictional force(s) Contact Point 1
    Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '|<=|', 'star', Opt);
    clear Bleft cLeft Bright cRight
    %friction cone clamping of frictional force(s) Contact Point 2
    [numConstraints, uIndex, Bleft, cLeft, Bright, cRight] = InitBC(1,[2, 4],[0 1],[],[obj.nu_p 0],[]); % friction cone clamping of frictional force(s) Contact Point 2
    Opt = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, lv1, '|<=|', 'star', Opt);
    clear Bleft cLeft Bright cRight
end


%Coordinate transform
function [xc] = coordinateTransformSC(obj, xs)
    %Extract variables from xs
    ribi = xs(1:2);
    theta = xs(3);
    ripi = xs(4:5);
    %Direction Cosine Matrices 
    Cbi = Helper.C3_2d(theta); 
    %Convert state to body frame
    rbbi = Cbi*ribi;
    rbpi = Cbi*ripi;
    %Build xc
    rbpb = rbpi - rbbi;
    ry = rbpb(2);
    %Output
    xc = [ribi;theta;ry];
end
%Coordinate transform
function [xs] = coordinateTransformCS(obj, xc)
    %Direction cosine matrices
    Cbi = Helper.C3_2d(xc(3));
    %Convert state to body frame
    rbpb = [-obj.a/2;xc(4)];
    ripb = Cbi'*rbpb;
    ripi = xc(1:2)+ripb;
    %Forces
    xs = [xc(1:2);xc(3);ripi;xc(3)];
end

function [xc,uc] = fullTransformSC(obj,xs,us) % xs = [x;y;theta;xp;yp;thetap] us = 3x1 pusher velocity (in world frame)
    s = Simulator();
    [dxs, Fgen] = s.lineSimulator(xs,us);
    xc = obj.coordinateTransformSC(xs); % xc = [x;y;theta;ry] 
    uc = [Fgen(1:2); Fgen(3)-Fgen(4); Fgen(5)-Fgen(6); Fgen(7)*sign(Fgen(3)-Fgen(4))]; % uc = [fn1,fn2,ft1,ft2,dry]
end
function [xs,us] = fullTransformCS(obj, xc,uc) %   xc = [x;y;theta;ry] uc = [fn1,fn2,ft1,ft2,dry]
    us = obj.force2Velocity(xc,uc); % us = 3x1 pusher velocity (in world frame)
    xs = obj.coordinateTransformCS(xc); % xs = [x;y;theta;xp;yp;thetap]. (pose of slider AND pose of pusher in inertial coordinates)
end

function [dxs, F_gen] = lineSimulator(obj,xs,us)
    % Quasi-Static Pushing Simulator 
    %
    % USAGE
    %   dxs = s.lineSimulator(xs,us);
    %
    %   Parameters:
    %   xs        -   6x1 pose matrix [x;y;theta;xp;yp;thetap]. (pose of slider AND pose of pusher in inertial coordinates)
    %   us        -   3x1 pusher velocity (in world frame)

    %Extract variables from xs
    ribi = xs(1:2);
    theta = xs(3);
    ripi = xs(4:5);
    thetap = xs(6);
    %Extract variables from us
    vipi = us(1:2);
    dthetap = us(3);
    %Direction Cosine Matrices 
    Cbi = Helper.C3_2d(theta);
    Cpi = Helper.C3_2d(thetap);   
    %Convert state to body frame
    rbbi = Cbi*ribi;
    rbpi = Cbi*ripi;
    vbpi = Cbi*vipi;
    %Find rx, ry: In body frame   
    rpap = [0;obj.lp/2];
    rpcp = [0;-obj.lp/2];
    rbap = Cbi*Cpi'*rpap;
    rbcp = Cbi*Cpi'*rpcp;
    rbab = rbpi + rbap - rbbi;
    rbcb = rbpi + rbcp - rbbi;
    rbai = rbpi + rbap;
    rbci = rbpi + rbcp;
    rax = rbab(1);
    ray = rbab(2);
    rcx = rbcb(1);
    rcy = rbcb(2);
    rb{1} = rbab;
    rb{2} = rbcb;
    %kinematics
    npa_b = [1;0];
    tpa_b = [0;1];
    Dpa_b = [tpa_b -tpa_b];
    Jpa = [1 0 -rbab(2);...
           0 1 rbab(1)];
    vpa = vbpi + Helper.cross3d(dthetap, rbap); 
    %kinematics
    npc_b = [1;0];
    tpc_b = [0;1];
    Dpc_b = [tpc_b -tpc_b];
    Jpc = [1 0 -rbcb(2);...
           0 1 rbcb(1)];
    vpc = vbpi + Helper.cross3d(dthetap, rbcp);
    %find distances
    rbzi = rbbi-[obj.a/2;0];
    rbaz = rbai-rbzi;
    rbcz = rbci-rbzi;
    d{1} = -rbaz'*npa_b;
    d{2} = -rbcz'*npc_b;
    %build useful matrices pt a
    N{1} = npa_b'*Jpa;
    L{1} = Dpa_b'*Jpa;
    E{1} = [1;1];
    aMat{1} = -npa_b'*vpa;
    bMat{1} = -Dpa_b'*vpa;
    %build useful matrices pt a
    N{2} = npc_b'*Jpc;
    L{2} = Dpc_b'*Jpc;
    E{2} = [1;1];
    aMat{2} = -npc_b'*vpc;
    bMat{2} = -Dpc_b'*vpc;
    %collision check
    counter = 1;
    index_vec = [];
    for lv1=1:2
        if d{lv1}<0.0001 && abs(rb{lv1}(2))<obj.a/2
            index_vec(counter) = lv1;
            counter = counter+1;
        end
    end
    %number if contact points
    m=length(index_vec);
    %concatenate matrices
    N_tot= [];
    L_tot=[];
    E_tot=[];
    a_tot = [];
    b_tot = [];
    nu_tot = [];
    F_gen = [];
    for lv1=1:m
        N_tot = [N_tot;N{index_vec(lv1)}];
        L_tot = [L_tot;L{index_vec(lv1)}];
        E_tot = blkdiag(E_tot,E{index_vec(lv1)});
        a_tot = [a_tot;aMat{index_vec(lv1)}];
        b_tot = [b_tot;bMat{index_vec(lv1)}];
        nu_tot = blkdiag(nu_tot,obj.nu_p);
    end
    %solve LCP
    if m
        %LCP matrices
        M = [N_tot*obj.A_ls*N_tot' N_tot*obj.A_ls*L_tot' zeros(m,m);...
            L_tot*obj.A_ls*N_tot' L_tot*obj.A_ls*L_tot' E_tot;...
            nu_tot -E_tot' zeros(m,m)];
        q = [a_tot;b_tot;zeros(m,1)];
        F_gen = LCP(M,q);
        fn = F_gen(1:m);
        ft = F_gen(m+1:m+1+2*m-1);
        twist = obj.A_ls*(N_tot'*fn + L_tot'*ft);
    else
        twist=[0;0;0];
    end
    %Compute twist (in body frame)
    twist_body_b = twist;
    %convert to inertial frame
    S_ib = [Cbi' [0;0];0 0 1];
    twist_body_i = S_ib*twist_body_b;
    %Return combined velocites of slider and pusher
    dxs = [twist_body_i; us];
end

end
    
end

classdef Controller < dynamicprops & PusherSliderSystem & Friction
    properties (Constant)
        h = 0.03;
        steps=50;
        num_xcStates = 4;
        num_ucStates = 5;
        uc_eq= [0.1634;0.1634;0;0;0]; %For Linearization purposes
        xc_eq = [0;0;0;0]; %For Linearization purposes
        Q_MPC = 10*diag([1,3,1,0]);
        Q_MPC_final= 2000*diag([1,3,1,0]);
        R_MPC = 1*diag([1,1,1,1,0]);
    end
    
    properties
        A_linear;
        B_linear;
        A_bar;
        B_bar;
        t_star;
        xc_star;
        uc_star
        xs_star;
        us_star;
    end
   
end

