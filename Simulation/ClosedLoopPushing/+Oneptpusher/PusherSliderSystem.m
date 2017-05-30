classdef PusherSliderSystem < dynamicprops 
    % Contains all object properties
    properties (Constant)
        %Pusher constants
        a=0.09;
        b=0.09;
        rho = 10000;
        height = 0.013;
    end
    
    properties
       m;
       V; 
       A; 
       d;
    end
   
    methods
        %% Constructor
        function obj = PusherSliderSystem()  
            %Set constant equations
            obj.A = obj.a*obj.b;
            obj.V = obj.A*obj.height;
            obj.m = obj.rho*obj.V;
        end
    end
end