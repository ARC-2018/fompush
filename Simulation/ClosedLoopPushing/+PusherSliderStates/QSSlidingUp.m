classdef QSSlidingUp < HybridInterfaces.HybridState
properties
    motion_function;
    controller_matrix;
    name = 'Sliding Up';
    x;
    u;
    E;
    D;
    g;
end
properties (Access = private)
    epsilon = 0.005;
    gamma_top;
    gamma_top_simplified;
    C_top;
    rx;
    rx2;
end
methods
    function obj = QSSlidingUp(a, nu_p, c2)
        obj.rx = a / 2.0;
        obj.rx2 = obj.rx * obj.rx;
        obj.x = sym('x', [1, 4]); % x, y, theta, ry
        obj.u = sym('u', [1, 2]); % v_n, v_t
        obj.gamma_top = (nu_p * c2 - obj.rx * obj.x(4) + nu_p * obj.rx2) / (c2 + obj.x(4)^2 - nu_p * obj.rx * obj.x(4));
        obj.gamma_top_simplified = matlabFunction(obj.gamma_top);
        Cbi = Helper.C3_2d(obj.x(3));
        Q = [c2 + obj.rx2 obj.rx * obj.x(4);
             obj.rx * obj.x(4) c2 + obj.x(4)^2] ./ (c2 + obj.rx2 + obj.x(4)^2);
        P = [1 0; obj.gamma_top 0];
        b = [(-obj.x(4) + obj.gamma_top * obj.rx) / (c2 + obj.rx2 + obj.x(4)^2) 0];
        c = [-obj.gamma_top 1];
        %Build nonlinear function
        obj.controller_matrix = [Cbi.' * Q * P; b; c];
        obj.motion_function = obj.controller_matrix * obj.u.';
        obj.C_top = jacobian(obj.gamma_top, obj.x);
        obj = obj.SetLinearConstraints();
        obj = obj.SetLinearMatrices();
    end
%     function [E, D, g] = GetLinearConstraints(obj, x_star, u_star)
%         C_top_subs = double(subs(obj.C_top, obj.x(4), x_star(4)));
%         gamma_top_star = double(subs(obj.gamma_top, obj.x(4), x_star(4)));
%         E = obj.u(1) * C_top_subs;
%         D = [gamma_top_star -1];
%         g = [u_star(2) - gamma_top_star * u_star(1) - obj.epsilon];
%     end
    
    function obj = SetLinearConstraints(obj)
        obj.E = matlabFunction(obj.u(1) * obj.C_top);
        obj.D = matlabFunction([obj.gamma_top -1]);
        obj.g = matlabFunction([obj.u(2) - obj.gamma_top * obj.u(1) - obj.epsilon]);
    end
    
    function [A, B, F, D, E, g] = GetLinearMatrices(obj, x, u)
        A = obj.A(u(1), x(3), x(4));
        B = obj.B(x(3), x(4));
        F = obj.F(u(1), u(2), x(3), x(4));
        E = obj.E(u(1), x(4));
        D = obj.D(x(4));
        g = obj.g(u(1), u(2), x(4));
    end
    
    function [F] = GetMotionFunction(obj, x, u)
        F = obj.F(u(1), u(2), x(3), x(4));
    end
    
    function [B_state, D, g] = GetInitialStateMatrices(obj, x_0_state, u_0_star)
        B_state = obj.B(x_0_state(3), x_0_state(4));
        D = obj.D(x_0_state(4));
        g = obj.g(u_0_star(1), u_0_star(2), x_0_state(4));
    end
    
    function is_inside = CheckConstraints(obj, x, u)
        gamma = u(2) / u(1);
        gamma_top_subs = obj.gamma_top_simplified(x(4));
        is_inside = gamma > gamma_top_subs;
    end
end 
end