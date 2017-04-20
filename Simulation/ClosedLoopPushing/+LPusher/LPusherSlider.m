classdef LPusherSlider < dynamicprops
    % Contains all object properties
    properties (Constant)
        %Pusher constants
        a=0.09;
        b=0.09;
        lp=0.030;
        rho = 10000;
        height = 0.013;
        A = LPusher.LPusherSlider.a * LPusher.LPusherSlider.b;
        V = LPusher.LPusherSlider.A * LPusher.LPusherSlider.height;
        m = LPusher.LPusherSlider.rho * LPusher.LPusherSlider.V;
        d = LPusher.LPusherSlider.lp / 2;
    end
   
end