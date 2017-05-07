classdef QSPusherSlider
%QSPUSHERSLIDER Class which contains all the variables of the Quasi Static
%Pusher Slider model as well as methods to obtain its Animation Data
properties (Constant)
    %Pusher constants
    a = 0.09;
    b = 0.09;
    diag = sqrt(Models.QSPusherSlider.a^2 + Models.QSPusherSlider.b^2);
    A = Models.QSPusherSlider.a * Models.QSPusherSlider.b;
    V = Models.QSPusherSlider.A * Models.QSPusherSlider.height;
    nu = 0.35;
    nu_pusher = 0.3;
    rho = 10000;
    m = Models.QSPusherSlider.rho * Models.QSPusherSlider.V;
    height = 0.013;
    f_max = (Models.QSPusherSlider.nu * Models.QSPusherSlider.m * Helper.g);
    m_max = Models.QSPusherSlider.m_max_funct(Models.QSPusherSlider.nu, Models.QSPusherSlider.m);
    c = Models.QSPusherSlider.m_max / Models.QSPusherSlider.f_max; 
    m_pert = Models.QSPusherSlider.m * 1.0; 
    nu_pert = Models.QSPusherSlider.nu * 1.0; 
    nu_pusher_pert = Models.QSPusherSlider.nu_pusher * 1.0; 
    f_max_pert = (Models.QSPusherSlider.nu_pert * Models.QSPusherSlider.m_pert * Helper.g);
    m_max_pert = Models.QSPusherSlider.m_max_funct(Models.QSPusherSlider.nu_pert, Models.QSPusherSlider.m_pert);
    c_pert = Models.QSPusherSlider.m_max_pert / Models.QSPusherSlider.f_max_pert;
end
methods (Static, Access = private)
    
function n_f = m_max_funct(nu, m)     
    n_f_integrand = @(p1, p2) (nu * m * Helper.g / Models.QSPusherSlider.A) * sqrt([p1; p2; 0]' * [p1; p2; 0]);
    n_f = Helper.DoubleGaussQuad(n_f_integrand, -Models.QSPusherSlider.a / 2, Models.QSPusherSlider.a / 2, -Models.QSPusherSlider.b / 2, Models.QSPusherSlider.b / 2);
end     
end
%% Animation Data
methods (Static)
    
function [x_slider, y_slider, x_pusher, y_pusher] = GetPusherSliderPolygons(x_state)
    %Declare variables
    d = x_state(4);
    theta = x_state(3);
    %Find distance d
    Rbi = Helper.C3_2d(x_state(3));
    ribi = x_state(1:2);
    R = 0.0045; % Radius of the Pusher Circle
    rbpb = [-(Models.QSPusherSlider.a/2 + R); d];
    ripb = Rbi.' * rbpb;
    ripi = ripb + ribi;
    %% Vertices of Slider Polygon
    P1 = x_state(1:2) + Rbi.' * [-Models.QSPusherSlider.a/2; -Models.QSPusherSlider.b/2];
    P2 = x_state(1:2) + Rbi.' * [Models.QSPusherSlider.a/2; -Models.QSPusherSlider.b/2];
    P3 = x_state(1:2) + Rbi.' * [Models.QSPusherSlider.a/2; Models.QSPusherSlider.b/2];
    P4 = x_state(1:2) + Rbi.' * [-Models.QSPusherSlider.a/2; Models.QSPusherSlider.b/2];
    %Build vector of vertices
    x_slider = [P1(1) P2(1) P3(1) P4(1)];
    y_slider = [P1(2) P2(2) P3(2) P4(2)];
    pusher_center = ripi;
    %% Points of Pusher Circle
    numPoints = 100; %Number of points making up the circle
    theta_vec = linspace(0, 2 * pi, numPoints);
    x_pusher = pusher_center(1) + R * cos(theta_vec);
    y_pusher = pusher_center(2) + R * sin(theta_vec);
end

function [x_lb, x_up, y_lb, y_ub] = GetPlotLimits(x_trajectory)
    x_lb = min(x_trajectory(1, :)) - Models.QSPusherSlider.diag * .6;
    x_up = max(x_trajectory(1, :)) + Models.QSPusherSlider.diag * .6;
    y_lb = min(x_trajectory(2, :)) - Models.QSPusherSlider.diag * .6;
    y_ub = max(x_trajectory(2, :)) + Models.QSPusherSlider.diag * .6;
end

end
end

