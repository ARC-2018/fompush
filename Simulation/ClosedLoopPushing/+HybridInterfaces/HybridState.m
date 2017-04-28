classdef HybridState < matlab.mixin.Heterogeneous & handle
% Interface that defines the properties and methods required for a hybrid
% state
properties (Abstract) % Abstract properties cannot define access methods,
   % they have to be defined in the implementation
   motion_function; % Kinematics function
   name;
   x;
   u;
   E;
   D;
   g;
end
properties
   A;
   B;
   F;
end
% A{lv1} = subs(A{lv1},{x,y,theta ry},{obj.x_eq(1),obj.x_eq(2),obj.x_eq(3), obj.ry_star})
methods
    function obj = SetLinearMatrices(obj)
        obj.F = matlabFunction(obj.motion_function);
        A_symbolic = jacobian(obj.motion_function, obj.x);
        obj.A = matlabFunction(A_symbolic);
        B_symbolic = jacobian(obj.motion_function, obj.u);
        obj.B = matlabFunction(B_symbolic);
        obj = obj.SetLinearConstraints();
    end
end
methods (Abstract)
    [E, D, g] = SetLinearConstraints(obj);
    [A, B, D, E, g] = GetLinearMatrices(obj, x_state, u_state);
    [F] = GetMotionFunction(obj, x_state, u_state);
    [B, D, g] = GetInitialStateMatrices(obj, x_state0, u_star0);
    % Checks whether the state constraints are fulfilled
    is_inside = CheckConstraints(obj, x_state, u_state);
end
    
end

