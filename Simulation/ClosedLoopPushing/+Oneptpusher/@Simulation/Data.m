%% Animation Data
function AniOutput = Data(obj, ite, Dataset)  
%         Dataset    
    x_state = obj.xs{Dataset};
    u_state = obj.us{Dataset};
    %Declare variables
    x     = x_state(ite, 1);
    y     = x_state(ite, 2);
    theta = x_state(ite, 3);
    xp     = x_state(ite, 4);
    yp     = x_state(ite, 5);
    u1 = u_state(ite, 1);
    u2 = u_state(ite, 2);
    %Rotation matrix
    Rb = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    %kinematics
    Cbi = Helper.C3_2d(theta);
    %position vectors (inertial frame)
    ribi = [x;y];
    ripi = [xp;yp];
    ripb = ripi-ribi;
    %position vectors (other frames)
    rbpb = Cbi*ripb;
    rbbi = Cbi*ribi;
    rbpi = Cbi*ripi;
    %find distances
    rbzi = rbbi-obj.a/2;
    rbpz = rbpi-rbzi;
    n_b = [1;0];
    
    rb{1} = rbpb;
    
    %% Calculate distances
    d{1} = -rbpz'*n_b;
    %% contact boolean variables
    contact = {};
    for lv1=1:1
        if d{lv1}<0.001 && abs(rb{lv1}(2))<obj.a/2
            contact{lv1} = 1;
        else
            contact{lv1} = 0;
        end
    end
%     [d{1} d{2};contact{1} contact{2}]
    %% Slider square
    % Vertices of Polygon (slider)
    P1b = [x; y] + Rb*[-obj.a/2; -obj.b/2];
    P2b = [x; y] + Rb*[obj.a/2; -obj.b/2];
    P3b = [x; y] + Rb*[obj.a/2; obj.b/2];
    P4b = [x; y] + Rb*[-obj.a/2; obj.b/2];
    %Build vector of vertices (slider)
    x1b = [P1b(1) P2b(1) P3b(1) P4b(1)];
    y1b = [P1b(2) P2b(2) P3b(2) P4b(2)];
    %% Pusher circles geometry
    Radius = 0.0035;
    numPoints=100;
    theta_vec=linspace(0,2*pi,numPoints); %100 evenly spaced points between 0 and 2pi
    rho=ones(1,numPoints)*Radius; %Radius should be 1 for all 100 points
    [X,Y] = pol2cart(theta_vec,rho); 
    %Point p (pusher)
    posp = ripi + Rb*([-Radius;0]);
    X_circle_p = X+posp(1);
    Y_circle_p = Y+posp(2);
    %Set ouputs
    AniOutput.x1b = x1b; 
    AniOutput.y1b = y1b; 
    AniOutput.X_circle_p = X_circle_p;
    AniOutput.Y_circle_p = Y_circle_p;
    AniOutput.posp = posp;
    AniOutput.contact = contact;
end