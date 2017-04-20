classdef LPController < LPusher.Friction & LPusher.LPusherSlider & matlab.mixin.Heterogeneous

properties (Constant)
    num_x = 4;
    num_u = 5;
    uc_eq= [0.1634;0.1634;0;0;0]; %For Linearization purposes
    xc_eq = [0;0;0;0]; %For Linearization purposes
end

properties
    
    A_linear;
    B_linear;
    t_star;
    xc_star;
    uc_star
    xs_star;
    us_star;
end

methods
    function obj = LPController()
        obj.symbolicLinearize;
    end
    
    function [Opt, H, Ai, bi, Ae, be] = buildConstraints(obj, Opt, index, x0, x_c, u_c, num_steps, h, Q_MPC, R_MPC, H)
        %Add dynamic constraints
        [Ai, bi, Ae, be, H] = obj.buidDynConstraintsAndCost(Opt, index, x0 - x_c, num_steps, h, Q_MPC, R_MPC, H);
        %Add mode independent constraints
        [A1,b1,A2,b2] = obj.buidModeIndep(Opt, index);
        Ai = [Ai; A1];
        bi = [bi; b1];
        Ae = [Ae; A2];
        be = [be; b2];
        [A1,b1,A2,b2] = obj.buildModeDepConstraints(Opt, index);
        Ai = [Ai; A1];
        bi = [bi; b1];
        Ae = [Ae; A2];
        be = [be; b2];
    end
    
    function [Ai, bi, Ae, be, H] = buidDynConstraintsAndCost(obj, Opt, index, delta_x, num_steps, h, Q_MPC, R_MPC, H)
        Ai = [];
        bi = [];
        Ae = [];
        be = [];
        numConstraints = obj.num_x;
        Aeq = zeros(numConstraints, Opt.nv);
        %Special case of initial conditions
        A_bar = eye(4) + h * obj.A_linear;
        beq = zeros(obj.num_x,1);
        if index ~=1
            Aeq(:,Opt.vars.x.i(1:obj.num_x,index-1))= -A_bar;
        end
        Aeq(:,Opt.vars.x.i(1:obj.num_x,index))  = eye(obj.num_x);
        B_bar = h * obj.B_linear;
        Aeq(:,Opt.vars.u.i(1:obj.num_u,index))=  -B_bar;
        if index == 1 % it has to be here to not add lines if unnecessary  
            beq(1:4) = [A_bar*(delta_x)];
        end
        Ae = [Ae; Aeq];
        be = [be; beq];
        H(Opt.vars.x.i(1:length(Q_MPC), index), Opt.vars.x.i(1:length(Q_MPC), index)) = Q_MPC;
        H(Opt.vars.u.i(1:length(R_MPC), index), Opt.vars.u.i(1:length(R_MPC), index)) = R_MPC;
        %Final Cost
        if index == num_steps
            Q_final = 10 * dare(A_bar, B_bar, Q_MPC, R_MPC);
            H(Opt.vars.x.i(1:obj.num_x, index), Opt.vars.x.i(1:obj.num_x, index)) = Q_final;
        end
    end
    
    function [Ai, bi, Ae, be] = buidModeIndep(obj, Opt, index)
        Ai = [];
        bi = [];
        Ae = [];
        be = [];
        %positive normal force(s)
        numConstraints = 2;
        uIndex = [1,2];
        Bleft = [];
        cLeft  = [];
        Bright  = eye(2);
        cRight = [];
        [A1,b1,A2,b2] = obj.addPusherConstraints2(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '<=', 'star', Opt);
        Ai = [Ai; A1];
        bi = [bi; b1];
        Ae = [Ae; A2];
        be = [be; b2];
        %friction cone clamping of frictional force(s) Contact Point 1
        numConstraints = 1;
        uIndex = [1,3];
        Bleft = [0 1];
        cLeft  = [];
        Bright  = [obj.nu_p 0];
        cRight = [];
        [A1,b1,A2,b2] = obj.addPusherConstraints2(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '|<=|', 'star', Opt);
        Ai = [Ai; A1];
        bi = [bi; b1];
        Ae = [Ae; A2];
        be = [be; b2];
        %friction cone clamping of frictional force(s) Contact Point 2
        numConstraints = 1;
        uIndex = [2,4];
        Bleft = [0 1];
        cLeft  = [];
        Bright  = [obj.nu_p 0];
        cRight = [];
        [A1,b1,A2,b2] = obj.addPusherConstraints2(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '|<=|', 'star', Opt); 
        Ai = [Ai; A1];
        bi = [bi; b1];
        Ae = [Ae; A2];
        be = [be; b2];
    end
    
    % Helper function Frank uses
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
        % TODO: CHANGE
        for index=1:2 %index represents contact point 1 and 2
            rb{index} = [rx*1;ry*1]+sign_vec(index)*[0;obj.d*1]; %position of contact point index
            Jb{index} = [1 0 -rb{index}(2);... 
                         0 1 rb{index}(1)];
            for lv2=1:2 %lv2 represents left of right border (FC/MC)
                n{index} = [1;0];
                t{index} = [0;1];
                N{index} = transpose(n{index})*Jb{index};
                T{index} = transpose(t{index})*Jb{index};
            end
        end
        %Motion equations (nonlinear)
        N_tot = [N{1};N{2}];
        T_tot = [T{1};T{2}];
        Vb = obj.A_ls*(transpose(N_tot) * fn + transpose(T_tot)*ft);
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
    end
    
    function [Ai, bi, Ae, be] = addPusherConstraints2(obj, Bleft, Bright, cLeft, cRight, numConstraints, I, index, equalFlag, nominalFlag, Opt)
    %       if strcmp(flag,'leq')
        Ai = [];
        bi = [];
        Ae = [];
        be = [];
        if isempty(Bleft)
            Bleft = zeros(numConstraints);
        end
        if isempty(Bright)
            Bright = zeros(numConstraints);
        end
        if isempty(cLeft)
            cLeft = zeros(numConstraints,1);
        end
        if isempty(cRight)
            cRight = zeros(numConstraints,1);
        end
        %% If lower than equal condition
        if ~strcmp(equalFlag,'|<=|') || ~strcmp(equalFlag,'|>=|')
            %Compute A
            A = Bleft-Bright;
            b = cRight-cLeft;
            Amatrix = zeros(numConstraints, Opt.nv);
            Amatrix(:,Opt.vars.u.i(I,index)) = A;
            %Compute b
            if strcmp(nominalFlag,'star')
                u_star = [];
                for lv3=1:length(I)
                    u_star = [u_star;obj.uc_eq(I(lv3))];
                end
                bmatrix = - A*u_star;
            else
                bmatrix = b;
            end
            %Add constraint
            epsilon=0.0001;
            if strcmp(equalFlag,'<=')
                Ai = [Ai; Amatrix];
                bi = [bi; bmatrix];
            elseif strcmp(equalFlag,'<')
                Ai = [Ai; Amatrix];
                bi = [bi; bmatrix-epsilon];
            elseif strcmp(equalFlag,'>=')
                Ai = [Ai; -Amatrix];
                bi = [bi; -bmatrix];
            elseif strcmp(equalFlag,'>')
                Ai = [Ai; -Amatrix];
                bi = [bi; -bmatrix-epsilon]; %TODO: Change
            elseif strcmp(equalFlag,'==')
                Ae = [Ae; Amatrix];
                be = [be; bmatrix];
            end
        end
        %% If lower than absolute value condition
        if strcmp(equalFlag,'|<=|')
            aSignVec = [-1,1];
            eqSignVec = [1,-1];

            for lv4=1:2
                %Compute A
                A = eqSignVec(lv4)*(Bleft+aSignVec(lv4)*Bright);
                b = (cRight+aSignVec(lv4)*cLeft);
                Ain = zeros(numConstraints, Opt.nv);
                Ain(:,Opt.vars.u.i(I,index)) = A;
                %Compute b
                if strcmp(nominalFlag,'star')
                    u_star = [];
                    for lv3=1:length(I)
                        u_star = [u_star;obj.uc_eq(I(lv3))];
                    end
                    bin = -A*u_star;

                else
                    bin = b;
                    disp(bin);
                end
                Ai = [Ai; Ain];
                bi = [bi; bin];
            end
        end
    end

    function Opt = addPusherConstraints(obj, Bleft, Bright, cLeft, cRight, numConstraints, I, index, equalFlag, nominalFlag, Opt)
    %       if strcmp(flag,'leq')
        if isempty(Bleft)
            Bleft = zeros(numConstraints);
        end
        if isempty(Bright)
            Bright = zeros(numConstraints);
        end
        if isempty(cLeft)
            cLeft = zeros(numConstraints,1);
        end
        if isempty(cRight)
            cRight = zeros(numConstraints,1);
        end
        %% If lower than equal condition
        if ~strcmp(equalFlag,'|<=|') || ~strcmp(equalFlag,'|>=|')
            %Compute A
            A = Bleft-Bright;
            b = cRight-cLeft;
            Amatrix = zeros(numConstraints, Opt.nv);
            Amatrix(:,Opt.vars.u.i(I,index)) = A;
            %Compute b
            if strcmp(nominalFlag,'star')
                u_star = [];
                for lv3=1:length(I)
                    u_star = [u_star;obj.uc_eq(I(lv3))];
                end
                bmatrix = - A*u_star;
            else
                bmatrix = b;
            end
            %Add constraint
            epsilon=0.0001;
            if strcmp(equalFlag,'<=')
                Opt =  Opt.addLinearConstraints(Amatrix, bmatrix, [], []);
            elseif strcmp(equalFlag,'<')
                Opt =  Opt.addLinearConstraints(Amatrix, bmatrix-epsilon, [], []);
            elseif strcmp(equalFlag,'>=')
                Opt =  Opt.addLinearConstraints(-Amatrix, -bmatrix, [], []);
            elseif strcmp(equalFlag,'>')
                Opt =  Opt.addLinearConstraints(-Amatrix, -bmatrix-epsilon, [], []);
            elseif strcmp(equalFlag,'==')
                Opt =  Opt.addLinearConstraints([], [], Amatrix, bmatrix);
            end
        end
        %% If lower than absolute value condition
        if strcmp(equalFlag,'|<=|')
            aSignVec = [-1,1];
            eqSignVec = [1,-1];

            for lv4=1:2
                %Compute A
                A = eqSignVec(lv4)*(Bleft+aSignVec(lv4)*Bright);
                b = (cRight+aSignVec(lv4)*cLeft);
                Ain = zeros(numConstraints, Opt.nv);
                Ain(:,Opt.vars.u.i(I,index)) = A;
                %Compute b
                if strcmp(nominalFlag,'star')
                    u_star = [];
                    for lv3=1:length(I)
                        u_star = [u_star;obj.uc_eq(I(lv3))];
                    end
                    bin = -A*u_star;

                else
                    bin = b;
                    disp(bin);
                end
                Opt =  Opt.addLinearConstraints(Ain, bin, [], []);
            end
        end
    end
    
    function [xc, uc] = Simulator2Controller(obj, xs, us)
        % Transform coordinates from simulator state coordinates to controller
        %Extract variables from xs
        ribi = xs(1:2);
        theta = xs(3);
        ripi = xs(4:5);
        thetap = xs(6);
        %Direction Cosine Matrices 
        Cbi = Helper.C3_2d(theta);
        Cpi = Helper.C3_2d(thetap);   
        %Convert state to body frame
        rbbi = Cbi*ribi;
        rbpi = Cbi*ripi;
        %Build xc
        rbpb = rbpi - rbbi;
        ry = rbpb(2);
        %Output
        xc = [ribi;theta;ry];
        if nargin > 1
            uc = [0.1634;0.1634;0;0;0]; %TODO: Frank should solve this
        else
            uc = nan;
        end
    end
    
end
    
end

