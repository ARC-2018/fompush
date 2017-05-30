classdef Controller < dynamicprops & Oneptpusher.PusherSliderSystem & Oneptpusher.Friction
    properties (Constant)
        h = 0.03;
        %Optimization Program
        h_opt = 0.03;
        NumFam = 3;
%         steps = 2;
        steps=50;
        num_xcStates = 4;
        num_ucStates = 3;
        num_xsStates = 5;
        num_usStates = 2;
        uc_eq= [0.3268;0;0]; %For Linearization purposes
        xc_eq = [0;0;0;0]; %For Linearization purposes
        Q = 10*diag([1,3,1,0]);
        Qf= 2000*diag([1,3,1,0]);
        R = 1*diag([1,1,0]);
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
        Opt;
        A_fun;
        B_fun;
        s;
    end
   
    methods
        %% Constructor
        function obj = Controller()
            obj.s = Oneptpusher.Simulator();
            obj.symbolicLinearize;
%             obj.buildNominalTrajectory;
%             obj.buildProgram();
        end
        %Coordinate transform
        function [xc] = coordinateTransformSC(obj, xs)
            % Transform coordinates from simulator state coordinates to controller
            %  state coordinates 
            %
            %   USAGE
            %   [xc] = coordinateTransform(xs)
            %
            %   Parameters:
            %   xs        -   5x1 simulator pose matrix [x;y;theta;xp;yp]. (pose of slider AND pose of pusher in inertial coordinates)
            %   xc        -   4x1 controller state [x;y;theta;ry]. (pose of slider AND relative pose of center of pusher)
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
            % Transform coordinates from controller state coordinates to
            % simulator states coordinates
            %
            %   USAGE
            %   [xs] = coordinateTransform(xc)
            %
            %   Parameters:
            %   xs        -   5x1 simulator pose matrix [x;y;theta;xp;yp]. (pose of slider AND pose of pusher in inertial coordinates)
            %   xc        -   4x1 controller state [x;y;theta;ry]. (pose of slider AND relative pose of center of pusher)
            %Extract variables
            ribi = xc(1:2);
            theta = xc(3);
            ry = xc(4);
            %Direction cosine matrices
            Cbi = Helper.C3_2d(theta);
            %Convert state to body frame
            rbbi = Cbi*ribi;
            rbpb = [-obj.a/2;ry];
            ripb = Cbi'*rbpb;
            ripi = ribi+ripb;
            %Forces
            xs = [ribi;theta;ripi];
        end
        
        function [xc,uc] = fullTransformSC(obj,xs,us); 
            % Transform coordinates and inputs states from simulator to
            % controller description
            %  state coordinates 
            %
            %   USAGE
            %   [xc,uc] = fullTransformSC(xs,us)
            %
            %   Parameters:
            %   xs        -   5x1 simulator pose matrix [x;y;theta;xp;yp]. (pose of slider AND pose of pusher in inertial coordinates)
            %   xc        -   4x1 controller state [x;y;theta;ry]. (pose of slider AND relative pose of center of pusher)
            %   us        -   2x1 pusher velocity (in world frame)
            %   uc        -   3x1 input matrix [fn1,ft1,dry]
            [dxs, Fgen] = obj.s.pointSimulator(xs,us);
            xc = obj.coordinateTransformSC(xs);
            uc = [Fgen(1); Fgen(2)-Fgen(3); Fgen(4)*sign(Fgen(2)-Fgen(3))];
        end
        
        function [xs,us] = fullTransformCS(obj, xc,uc); 
            % Transform coordinates and inputs states from controller to
            % simulator description
            %
            %   USAGE
            %   [xs,us] = fullTransformCS(xc,uc)
            %
            %   Parameters:
            %   xs        -   5x1 simulator pose matrix [x;y;theta;xp;yp]. (pose of slider AND pose of pusher in inertial coordinates)
            %   xc        -   4x1 controller state [x;y;theta;ry]. (pose of slider AND relative pose of center of pusher)
            %   us        -   2x1 pusher velocity (in world frame)
            %   uc        -   3x1 input matrix [fn1,ft1,dry]
            us = obj.force2Velocity(xc,uc);
            xs = obj.coordinateTransformCS(xc);
        end
        
       
        %
        function [us] = force2Velocity(obj, xc, uc)
            % Transform coordinates and inputs from controller state coordinates to
            % simulator states coordinates
            %
            %   USAGE
            %   [us] = force2Velocity(obj, xc, uc)
            %
            %   Parameters:
            %   us        -   2x1 pusher velocity (in world frame)
            %   uc        -   3x1 input matrix [fn1,ft1,dry]
            %   xc        -   4x1 controller state [x;y;theta;ry]. (pose of slider AND relative pose of center of pusher)
            %Extract variables
            ribi = xc(1:2);
            theta = xc(3);
            ry = xc(4);
            %Direction cosine matrices
            Cbi = Helper.C3_2d(theta);
            %Convert state to body frame
            rbbi = Cbi*ribi;
            rbpb = [-obj.a/2;ry];
            ripb = Cbi'*rbpb;
            %Forces
            fn = uc(1);
            ft = uc(2);
            dry = uc(3);
            drbpb = [0;dry];
            %kinematics
            npa_b = [1;0];
            tpa_b = [0;1];
            Tpa_b = [tpa_b];
            Jpa = [1 0 -rbpb(2);...
                   0 1 rbpb(1)];
            %build useful matrices pt a
            N{1} = npa_b'*Jpa;
            T{1} = Tpa_b'*Jpa;
            %concatenate matrices
            N_tot= [];
            T_tot=[];
            for lv1=1:1
                N_tot = [N_tot;N{lv1}];
                T_tot = [T_tot;T{lv1}];
            end
            %Compute twist (in body frame)
            twist_body_b = obj.A_ls*(N_tot'*fn+T_tot'*ft);
            %convert to inertial frame
            S_ib = [Cbi' [0;0];0 0 1];
            twist_body_i = S_ib*twist_body_b;
            %Compute twist of pusher
            vel_pusher_i = zeros(2,1);
            vel_pusher_i(1:2) = twist_body_i(1:2) + Cbi'*drbpb + Helper.cross3d(twist_body_i(3), ripb);
            %output
            us = vel_pusher_i;
        end
        %%
        function [A, B] = motionConstraintMatrices(obj, xn, un)
            % Returns dynamic constraint matrices
            % condition and time index
            %
            %   USAGE
            %   [A, B] = dynamicConstraintMatrices(obj, xn, un)
            %   Form of equation: ddelta_xc = A * delta_xc + B*delta*uc
            %
            %   Parameters:
            %   xn        -   4x1 controller state [x;y;theta;ry] of nominal trajectory. 
            %   un        -   3x1 controller inputs [fn1,ft1,dry] of nominal trajectory.
            A = obj.A_fun(xn, un);
            B = obj.B_fun(xn, un);
        end
        %solve MPC optimization problem
        function [Aeq, beq, Ain, bin] = forceConstraintMatrices(obj, ~, un, mode)
            % Returns dynamic constraint matrices
            % condition and time index
            %
            %   USAGE
            %   [Aeq, beq, Ain, bin] = forceConstraintMatrices(obj, xn, un, mode)
            %   Form of equations: Aeq*uc = beq, where uc = [fn1,ft1,dry]
            %   Form of equations: Ain*uc <= bin, where uc = [fn1,ft1,dry]
            %
            %   Parameters:
            %   xn        -   4x1 controller state [x;y;theta;ry] of nominal trajectory. 
            %   un        -   3x1 controller inputs [fn1,ft1,dry] of nominal trajectory.
            %   mode      -   0:stick, 1:slide up (left), 2:slide down (right).
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Initialize matrices
            Aeq = [];
            beq = [];
            Ain = [];
            bin = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Mode independent constraints %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% tangential velocity bounds
            num_constraints = 1;
            index = [3];
            Ain_tmp = 1;
            bin_tmp = 0.05; %upper limit (i.e. dry<0.05)
            Ain = [Ain;zeros(num_constraints, length(un))];
            bin = [bin;zeros(num_constraints, 1)];
            Ain(end-num_constraints+1:end,index) = Ain_tmp;
            bin(end-num_constraints+1:end) = bin_tmp;
            clear Ain_tmp bin_tmp
            %% tangential velocity bounds
            num_constraints = 1;
            index = [3];
            Ain_tmp = -1;
            bin_tmp = 0.05; %lower limit (i.e. dry>-0.05)
            Ain = [Ain;zeros(num_constraints, length(un))];
            bin = [bin;zeros(num_constraints, 1)];
            Ain(end-num_constraints+1:end,index) = Ain_tmp;
            bin(end-num_constraints+1:end) = bin_tmp;
            clear Ain_tmp bin_tmp
            % positive normal force
            num_constraints = 1;
            index = [1];
            Ain_tmp = -1;
%             un(index)
            bin_tmp = [0] - Ain_tmp*un(index);
            Ain = [Ain;zeros(num_constraints, length(un))];
%             zeros(num_constraints, 1)
            bin = [bin;zeros(num_constraints, 1)];
            Ain(end-num_constraints+1:end,index) = Ain_tmp;
            bin(end-num_constraints+1:end) = bin_tmp;
            clear Ain_tmp bin_tmp
            %%  normal force bounds
            %upper bound (applied on deviation only, aka only linear term)
            num_constraints = 1;
            index = [1];
            Ain_tmp = 1;
            max_normal = 0.05; %upper limit (i.e. delta_fn< 0.05)
            bin_tmp = [max_normal]; 
            Ain = [Ain;zeros(num_constraints, length(un))];
            bin = [bin;zeros(num_constraints, 1)];
            Ain(end-num_constraints+1:end,index) = Ain_tmp;
            bin(end-num_constraints+1:end) = bin_tmp;
            clear Ain_tmp bin_tmp
            %lower bound (applied on deviation only, aka only linear term)
            num_constraints = 1;
            index = [1];
            Ain_tmp = -1;
            max_normal = 0.1; %lower limit (i.e. delta_fn> -0.1)
            bin_tmp = [max_normal]; 
            Ain = [Ain;zeros(num_constraints, length(un))];
            bin = [bin;zeros(num_constraints, 1)];
            Ain(end-num_constraints+1:end,index) = Ain_tmp;
            bin(end-num_constraints+1:end) = bin_tmp;
            clear Ain_tmp bin_tmp
            %%  friction cone (contact point 1)
            %upper border
            num_constraints = 1;
            index = [1,2];
            Ain_tmp = [-obj.nu_p 1];
            bin_tmp =  - Ain_tmp*un(index);
            Ain = [Ain;zeros(num_constraints, length(un))];
            bin = [bin;zeros(num_constraints, 1)];
            Ain(end-num_constraints+1:end,index) = Ain_tmp;
            bin(end-num_constraints+1:end) = bin_tmp;
            clear Ain_tmp bin_tmp
            %lower border
            num_constraints = 1;
            index = [1,2];
            Ain_tmp = -[obj.nu_p 1];
            bin_tmp =  - Ain_tmp*un(index);
            Ain = [Ain;zeros(num_constraints, length(un))];
            bin = [bin;zeros(num_constraints, 1)];
            Ain(end-num_constraints+1:end,index) = Ain_tmp;
            bin(end-num_constraints+1:end) = bin_tmp;
            clear Ain_tmp bin_tmp
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Mode Dependent Constraints %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if mode==1
                %sticking constraints
                num_constraints = 1;
                index = [3];
                Aeq_tmp = 1;
                beq_tmp =  - Aeq_tmp*un(index);
                Aeq = [Aeq;zeros(num_constraints, length(un))];
                beq = [beq;zeros(num_constraints, 1)];
                Aeq(end-num_constraints+1:end,index) = Aeq_tmp;
                beq(end-num_constraints+1:end) = beq_tmp;
                clear Aeq_tmp beq_tmp
            elseif mode==2
                %sliding up
                num_constraints = 1;
                index = [3];
                Ain_tmp = -1;
                bin_tmp =  - Ain_tmp*un(index);
                Ain = [Ain;zeros(num_constraints, length(un))];
                bin = [bin;zeros(num_constraints, 1)];
                Ain(end-num_constraints+1:end,index) = Ain_tmp;
                bin(end-num_constraints+1:end) = bin_tmp;
                clear Ain_tmp bin_tmp
                %tangential force of contact point 1 must lie on bondary of FC
                num_constraints = 1;
                index = [1,2];
                Aeq_tmp = [-obj.nu_p 1];
                beq_tmp =  - Aeq_tmp*un(index);
                Aeq = [Aeq;zeros(num_constraints, length(un))];
                beq = [beq;zeros(num_constraints, 1)];
                Aeq(end-num_constraints+1:end,index) = Aeq_tmp;
                beq(end-num_constraints+1:end) = beq_tmp;
                clear Aeq_tmp beq_tmp
            elseif mode==3
                %sliding up
                num_constraints = 1;
                index = [3];
                Ain_tmp = 1;
                bin_tmp =  -Ain_tmp*un(index);
                Ain = [Ain;zeros(num_constraints, length(un))];
                bin = [bin;zeros(num_constraints, 1)];
                Ain(end-num_constraints+1:end,index) = Ain_tmp;
                bin(end-num_constraints+1:end) = bin_tmp;
                clear Ain_tmp bin_tmp
                %tangential force of contact point 1 must lie on bondary of FC
                num_constraints = 1;
                index = [1,2];
                Aeq_tmp = [obj.nu_p 1];
                beq_tmp =  - Aeq_tmp*un(index);
                Aeq = [Aeq;zeros(num_constraints, length(un))];
                beq = [beq;zeros(num_constraints, 1)];
                Aeq(end-num_constraints+1:end,index) = Aeq_tmp;
                beq(end-num_constraints+1:end) = beq_tmp;
                clear Aeq_tmp beq_tmp
            end
            
        end
        %Build cost matrix
        function Opt = buildCost(obj, Opt, lv1)
            %Initialize matrices
            H = zeros(Opt.nv, Opt.nv);
            % State Cost
            if lv1<obj.steps
                H(Opt.vars.x.i(1:length(obj.Q),lv1), Opt.vars.x.i(1:length(obj.Q),lv1)) = obj.Q;
            else
                H(Opt.vars.x.i(1:length(obj.Qf),lv1), Opt.vars.x.i(1:length(obj.Qf),lv1)) = obj.Qf;
            end
            %Inputs Cost
            H(Opt.vars.u.i(1:length(obj.R),lv1), Opt.vars.u.i(1:length(obj.R),lv1)) = obj.R;
            %Add cost to program
            Opt = Opt.addCost(H, [], []);
        end
        %Build dynamic constraints
        function Opt = addMotionConstraints(obj, Opt, lv1)
            A_nom = obj.A_fun(obj.xc_eq, obj.uc_eq);
            B_nom = obj.B_fun(obj.xc_eq, obj.uc_eq);
            A_bar = eye(obj.num_xcStates)+obj.h*A_nom;
            B_bar = obj.h*B_nom;

            numConstraints = obj.num_xcStates;
            Aeq = zeros(numConstraints, Opt.nv);
            %Special case of initial conditions
            if lv1 ~=1
            Aeq(:,Opt.vars.x.i(1:obj.num_xcStates,lv1-1))= -A_bar;
            end
            Aeq(:,Opt.vars.x.i(1:obj.num_xcStates,lv1))  = eye(obj.num_xcStates);
            Aeq(:,Opt.vars.u.i(1:obj.num_ucStates,lv1))=  -B_bar;
            beq = zeros(obj.num_xcStates,1);
            Opt = Opt.addLinearConstraints([], [], Aeq, beq);
        end
        %Build forces constraints
        function Opt = addForceConstraints(obj, Opt, lv1, mode)
            xn = obj.xc_eq;
            un = obj.uc_eq;
            
            [Aeq, beq, Ain, bin] = obj.forceConstraintMatrices(xn,un,mode);
%            %finite lateral velocity
            numConstraints = size(Ain,1);
            Amatrix = zeros(numConstraints, Opt.nv);
            bmatrix = zeros(numConstraints, 1);
            Amatrix(:,Opt.vars.u.i(:,lv1)) = Ain;
            bmatrix = bin;
            Opt =  Opt.addLinearConstraints(Amatrix, bmatrix, [], []);
            clear Amatrix bmatrix
            
            numConstraints = size(Aeq,1);
            Amatrix = zeros(numConstraints, Opt.nv);
            bmatrix = zeros(numConstraints, 1);
            Amatrix(:,Opt.vars.u.i(:,lv1)) = Aeq;
            bmatrix = beq;
            Opt =  Opt.addLinearConstraints([], [], Amatrix, bmatrix);
            clear Amatrix bmatrix   
        end
        %Build nominal trajectory
        function obj = buildNominalTrajectory(obj)
            %  Build 4 matrices that encode the nominal trajectory as a
            %  function of time
            %   Parameters:
            %   xs_star        -   Nx6 state matrix nominal trajectory in simulator coordinates
            %   us_star        -   Nx3 input matrix nominal trajectory in simulator coordinates
            %   xc_star        -   Nx4 state matrix nominal trajectory in controller coordinates
            %   uc_star        -   Nx5 input matrix nominal trajectory in controller coordinates
            v_eq = 0.05;
            tf = 20;
            t0 = 0;
            h_star = 0.01;
            N_star = (1/h_star)*(tf-t0);
            obj.t_star = zeros(N_star,1);
            for lv1=1:N_star
                %Define nominal values (simulator coordinates)
                xsStarTemp = [v_eq*obj.t_star(lv1) 0 0 v_eq*obj.t_star(lv1)-obj.a/2 0]';
                usStarTemp = [v_eq 0]';
                %Define nominal values (controller coordinates)
                xcStarTemp = obj.coordinateTransformSC(xsStarTemp);
                ucStarTemp = obj.uc_eq;
                %Build control matrices A and B (symbolic linearization of motion
                %equations)
                obj.xs_star(lv1,:) = xsStarTemp';
                obj.us_star(lv1,:) = usStarTemp';
                obj.xc_star(lv1,:) = xcStarTemp';
                obj.uc_star(lv1,:) = ucStarTemp';
                if lv1<N_star
                    obj.t_star(lv1+1)  = obj.t_star(lv1) + h_star;
                end
            end
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
            syms fx fy m fn1 ft1 ry_dot
            %% Build states       
            xc = [x;y;theta;ry];
            uc = [fn1;ft1; ry_dot];
            fn = [fn1];
            ft = [ft1];
            %% DCM Matrices
            Cbi = Helper.C3_2d(theta);
            %% Kinematics
            sign_vec = [1 -1]*1;
            rx = -obj.a/2;
            for lv1=1:1 %lv1 represents contact point 1 and 2
                rb{lv1} = [rx*1;ry*1]+sign_vec(lv1)*[0;0]; %position of contact point lv1
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
            N_tot = [N{1}];
            T_tot = [T{1}];
            Vb = obj.A_ls*(transpose(N_tot)*fn + transpose(T_tot)*ft );
            C_tilde = [Cbi 0;0 0 1];
            f_non1 = C_tilde*Vb;
            f_non2 = ry_dot;
            f_non = [f_non1;f_non2];
            %Linearization
            A = jacobian(f_non,xc);
            B = jacobian(f_non,uc);
            % Substitute equilibrium states
            obj.A_fun = matlabFunction(A, 'Vars', {[x; y; theta; ry], [fn1; ft1 ; ry_dot]});
            obj.B_fun = matlabFunction(B, 'Vars', {[x; y ;theta; ry], [fn1 ;ft1 ; ry_dot]});

        end
        %Get nominal trajectory values at time T
        function [xcStar, ucStar, xsStar, usStar] = getStateNominal(obj, t)
            % Get nominal trajectory values at time T
            %
            %   Parameters:
            %   t       -   Time at which the nominal state is evaluated
            vecDif = t - obj.t_star;
            [valDif, indexDif] = min(abs(vecDif));
            xcStar = obj.xc_star(indexDif,:)';
            ucStar = obj.uc_star(indexDif,:)'; 
            xsStar = obj.xs_star(indexDif,:)';
            usStar = obj.us_star(indexDif,:)'; 
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
            [xcStar, ucStar, usStar, usStar] = obj.getStateNominal(t);
            %Build error state vector
            delta_xc = [xc - xcStar];
            A_nom = obj.A_fun(obj.xc_eq, obj.uc_eq);
            B_nom = obj.B_fun(obj.xc_eq, obj.uc_eq);
            A_bar = eye(4)+obj.h*A_nom;
            B_bar = obj.h*B_nom;
            %Loop through family of modes
            fVal = [];
            for lv2=1:3
            %Build optimization program
%             obj.Opt{lv2} = obj.buildProgram();
                %Update initial conditions
                obj.Opt{lv2}.beq(1:4) = zeros(4,1);
                obj.Opt{lv2}.beq(1:4) = [A_bar*delta_xc];
                % Solve Opt Program   
                options = optimoptions('quadprog','Display','none');
                [obj.Opt{lv2}, solvertime{lv2}, fval{lv2}] = obj.Opt{lv2}.solve;
                out_delta_u{lv2} = obj.Opt{lv2}.vars.u.value';
                out_delta_x{lv2} = obj.Opt{lv2}.vars.x.value';
                fVal = [fVal; fval{lv2}];
            end
            %Find mode schedule with lowest cost
            [minFOM indexFOM] = min(fVal);
            %Return first element of control sequence
            delta_u = out_delta_u{indexFOM}(1,1:obj.num_ucStates)';
            %Add feedforward and feedback controls together
            uc = delta_u + ucStar;
        end
       %% Build optimization program and constraints
        function obj = buildProgram(obj)
            for lv1=1:3
                Opt = MixedIntegerConvexProgram(false);
                Opt = Opt.addVariable('x', 'C', [obj.num_xcStates, obj.steps], -1000*ones(obj.num_xcStates,obj.steps), 1000*ones(obj.num_xcStates,obj.steps));
                Opt = Opt.addVariable('u', 'C', [obj.num_ucStates, obj.steps], -1000*ones(obj.num_ucStates,obj.steps), 1000*ones(obj.num_ucStates,obj.steps));
                                %Loop through steps of MPC
                for lv3=1:obj.steps 
                    %Add cost
                    Opt = obj.buildCost(Opt, lv3);
                    %Add dynamic constraints
                    Opt = obj.addMotionConstraints(Opt, lv3);
                    %Add mode independent constraints
%                     Opt = obj.buildModeIndepConstraints(Opt, lv3);
                    if lv3==1
                        if lv1==1 
                            %Add mode dependant constraints
                            Opt = obj.addForceConstraints(Opt, lv3, 1);%1: Sticking, 2: Sliding left(up), 3: Slide right(down)
                        elseif lv1==2 
                            %Add mode dependant constraints
                            Opt = obj.addForceConstraints(Opt, lv3, 2);%1: Sticking, 2: Sliding left(up), 3: Slide right(down)
                        elseif lv1==3 
                            %Add mode dependant constraints
                            Opt = obj.addForceConstraints(Opt, lv3, 3);%1: Sticking, 2: Sliding left(up), 3: Slide right(down)
                        end
                    else
                        %Add mode dependant constraints
                        Opt = obj.addForceConstraints(Opt, lv3, 1);%1: Sticking, 2: Sliding left(up), 3: Slide right(down)
                    end  
                end
                obj.Opt{lv1} = Opt;
            end
        end
    end
end