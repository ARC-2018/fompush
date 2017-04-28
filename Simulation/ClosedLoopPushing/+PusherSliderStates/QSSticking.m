classdef QSSticking < HybridInterfaces.HybridState
properties % TODO: Look how to f*** make this private
    motion_function;
    controller_matrix; % TODO: UnitTest this is equal to B
    name = 'Quasi Static Sticking State';
    x;
    u;
    E;
    D;
    g;
end
properties
% properties (Access = private)
    gamma_top;
    gamma_top_simplified;
    gamma_bottom;
    gamma_bottom_simplified;
    C_top;
    C_bottom;
    rx;
    rx2;
end
methods
    function obj = QSSticking(a, nu_p, c2)
        obj.rx = a / 2.0;
        obj.rx2 = obj.rx * obj.rx;
        obj.x = sym('x', [4, 1]); % x, y, theta, ry
        obj.u = sym('u', [2, 1]); % v_n, v_t
        obj.gamma_top = (nu_p * c2 - obj.rx * obj.x(4) + nu_p * obj.rx2) / (c2 + obj.x(4)^2 - nu_p * obj.rx * obj.x(4));
        obj.gamma_top_simplified = matlabFunction(obj.gamma_top);
        obj.gamma_bottom = (-nu_p * c2 - obj.rx * obj.x(4) - nu_p * obj.rx2) / (c2 + obj.x(4)^2 + nu_p * obj.rx * obj.x(4));
        obj.gamma_bottom_simplified = matlabFunction(obj.gamma_bottom);
        Cbi = Helper.C3_2d(obj.x(3));
        Q = [c2 + obj.rx2 obj.rx * obj.x(4);
             obj.rx * obj.x(4) c2 + obj.x(4)^2] ./ (c2 + obj.rx2 + obj.x(4)^2);
        P = eye(2); % Not really needed, but keeps consistency with paper
        b = [-obj.x(4), obj.rx] / (c2 + obj.rx2 + obj.x(4)^2);
        c = [0 0];
        %Build nonlinear function
        obj.controller_matrix = [Cbi.' * Q * P; b; c];
        obj.motion_function = obj.controller_matrix * obj.u;
        obj.C_top = jacobian(obj.gamma_top, obj.x);
        obj.C_bottom = jacobian(obj.gamma_bottom, obj.x);
        obj = obj.SetLinearMatrices();
    end
    
    
    function obj = SetLinearConstraints(obj)
        obj.E = matlabFunction(obj.u(1) * [-obj.C_top; obj.C_bottom]);
        obj.D = matlabFunction([-obj.gamma_top 1; obj.gamma_bottom -1]);
        obj.g = matlabFunction([-obj.u(2) + obj.gamma_top * obj.u(1);
                                 obj.u(2) - obj.gamma_bottom * obj.u(1)]);
    end
%     function [E, D, g] = GetLinearConstraints(obj, x_star, u_star)
%         C_top_subs = double(subs(obj.C_top, obj.x(4), x_star(4)));
%         C_bottom_subs = double(subs(obj.C_bottom, obj.x(4), x_star(4)));
%         gamma_top_star = double(subs(obj.gamma_top, obj.x(4), x_star(4)));
%         gamma_bottom_star = double(subs(obj.gamma_bottom, obj.x(4), x_star(4)));
%         E = u_star(1) * [-C_top_subs; C_bottom_subs];
%         D = [-gamma_top_star 1; gamma_bottom_star -1];
%         g = [-u_star(2) + gamma_top_star * u_star(1);
%              u_star(2) - gamma_bottom_star * u_star(1)];
%     end
    
    function [A, B, D, E, g] = GetLinearMatrices(obj, x, u) % x = [x, y, theta, ry] u = [v_n, v_t]
        A = obj.A(u(1), u(2), x(3), x(4));
        B = obj.B(x(3), x(4));
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
        gamma_bottom_subs = obj.gamma_bottom_simplified(x(4));
        gamma_top_subs = obj.gamma_top_simplified(x(4));
        is_inside = gamma >= gamma_bottom_subs & gamma <= gamma_top_subs;
    end
end
end

