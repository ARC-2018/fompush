function [xc, uc] = Simulator2Controller(xs, us)
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