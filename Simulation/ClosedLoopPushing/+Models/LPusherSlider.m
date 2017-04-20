classdef LPusherSlider
    % Contains all object properties
    properties (Constant)
        %Pusher constants
        a=0.09;
        b=0.09;
        lp=0.030;
        rho = 10000;
        height = 0.013;
        A = Models.LPusherSlider.a * Models.LPusherSlider.b;
        V = Models.LPusherSlider.A * Models.LPusherSlider.height;
        m = Models.LPusherSlider.rho * Models.LPusherSlider.V;
        d = Models.LPusherSlider.lp/2;
        
        %Friction constants
        nu = 0.35;
        nu_p = 0.3;
    end
end