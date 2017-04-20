classdef Friction < dynamicprops
    % Class that defines the limit surface (A_ls), used in controller
    % design AND simulator
    properties (Constant)
        %Pusher constants
        nu = 0.35;
        nu_p = 0.3;
        f_max = (LPusher.Friction.nu * LPusher.LPusherSlider.m*Helper.g);
        m_max = LPusher.Friction.m_max_funct(LPusher.Friction.nu, LPusher.LPusherSlider.m);
        c = LPusher.Friction.m_max / LPusher.Friction.f_max; 
    end
    
    properties
       A_ls;
    end
    
    
    methods
        %% Constructor
        function obj = Friction()  
            %Compute f_max and m_max
            %Build limit surface representation
            syms fx fy m 
            F  = [fx;fy;m];
            %% 1) Limit surface
            H = fx^2/LPusher.Friction.f_max^2 + fy^2/LPusher.Friction.f_max^2 + m^2/LPusher.Friction.m_max^2;
            obj.A_ls = double(hessian(H,F));
        end
    end
    methods (Static)
        %Implement polar coordinates version
        function n_f = m_max_funct(nu, m)     
            n_f_integrand = @(p1,p2) (nu * m * Helper.g / LPusher.LPusherSlider.A) * sqrt([p1;p2;0]'*[p1;p2;0]);
            n_f = Helper.DoubleGaussQuad(n_f_integrand,-LPusher.LPusherSlider.a/2,LPusher.LPusherSlider.a/2,-LPusher.LPusherSlider.b/2, LPusher.LPusherSlider.b/2);
        end 
    end
end