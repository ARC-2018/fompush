classdef Animator
    %ANIMATOR Summary of this class goes here
    %   Detailed explanation goes here

properties (Constant)
    %Pusher constants
    a=0.09;
    b=0.09;
    d=0.015;
    nu = 0.35;
    nu_p = 0.3;
    rho = 10000;
    height = 0.013;
    %Euler Integration
    num_states = 6;
    %Optimization Program
    h_opt = 0.03;
    steps = 50; 
    NumFam = 3;
    num_int = 1;
    num_vars = 4;
    num_xvars = 3;
    num_dvars = 1;
    num_inputs = 5;
end    


properties
    t;
    xs;
    us;
    xc;
    uc;
    xs_state;
    us_state;
    xc_state;
    uc_state;
    uc_eq;
    xs_star_state;
    us_star_state;
    xc_star_state;
    uc_star_state;
    N;
    FileName;
    SimName;
    FilePath;
    NumSim;
    Ani;
end


methods

    %% Constructor
    function obj = Animator(SimName)  
        %% Save data in new folder
        obj.SimName = SimName;
        tempName = strcat('Test/',obj.SimName);
        mkdir(pwd,tempName);
        obj.FilePath = strcat(pwd,'/Test/',obj.SimName);
        obj.FileName = strcat(obj.FilePath,'/',obj.SimName);
        %Nominal Trajectory  #Frank hack
        %Initialize variables

    end
end

methods
function obj = plot(obj, flag)

    if flag == 0
        index = [1;2];
    elseif flag == 1
        index = [3,4];
    elseif flag == 2
        index = [5];
    end
        Name = 'test';
        Figures.(Name)=Figure;
        Figures.(Name).filename = Name;  
        Figures.(Name).Create(length(index),1); 

        if flag==0
            Figures.(Name).xData = {[obj.t, obj.t];[obj.t, obj.t]};
            Figures.(Name).yLabel={'$f_{n1}$ (N)';'$f_{n2}$ (N)'}; 
            Y1  = obj.uc{2}(:,index(1))*0;
            Y2  = obj.uc{2}(:,index(1))*0;
            Figures.(Name).yData = {[Y1,obj.uc{2}(:,index(1))];[Y2,obj.uc{2}(:,index(2))]};
            Figures.(Name).xLabel={'t(s)';'t(s)'};
            Figures.(Name).Color = {['k','b'];['k','r']};
        elseif flag==1
            Figures.(Name).xData = {[obj.t, obj.t, obj.t];[obj.t, obj.t, obj.t]};
            Figures.(Name).yLabel={'$f_{t1}$ (N)';'$f_{t2}$ (N)'}; 
            Y1  = obj.uc{2}(:,1)*obj.nu_p;
            Y2  = -obj.uc{2}(:,1)*obj.nu_p;
            Y3  = obj.uc{2}(:, 2)*obj.nu_p;
            Y4  = -obj.uc{2}(:,2)*obj.nu_p;
            Figures.(Name).yData = {[Y1,Y2,obj.uc{2}(:,index(1))];[Y3,Y4,obj.uc{2}(:,index(2))]};
            Figures.(Name).xLabel={'t(s)';'t(s)'};
            Figures.(Name).Color = {['k','k','b'];['k','k','r']};
        else
            Figures.(Name).xData = {[obj.t, obj.t]};
            Figures.(Name).yLabel={'$dot{r}_y$ (m/s)'}; 
            Y1  = obj.uc{2}(:,5)*0;
            Figures.(Name).yData = {[Y1,obj.uc{2}(:,index(1))]};
            Figures.(Name).xLabel={'t(s)'};
            Figures.(Name).Color = {['k','b']};
        end

        Figures.(Name).Title = {Name};               
        Figures.(Name).Plot2d;
    %     Figures.(Name).Save(Foldername);
end

%% Main Animation 
function obj = Animate(obj, flag)
%Initialize figure
obj.Ani = figure('Color', 'w', 'OuterPosition', [0, 0, 960, 1080], 'PaperPosition', [0, 0, 11, (6/8)*11]);
%figure properties
font_size  = 25;
line_size  = 15;
line_width = 2;
set(gcf,'Renderer','OpenGL');
set(gca,'FontSize',20)
axis equal
xlabel('x(m)','fontsize',font_size,'Interpreter','latex', 'FontSize', font_size);
ylabel('y(m)','fontsize',font_size,'Interpreter','latex', 'FontSize', font_size);
xlim([-.1 .4]);
ylim([-0.1 0.15]);
%Animation parameters
tf = obj.t(end);
N = length(obj.t);
accFactor = 5;
x_state = obj.xs{1};
%create movie file
videoname = strcat(obj.FilePath,'/',(obj.SimName),'.avi');
v = VideoWriter(videoname);
fps = int64(N/(accFactor*tf));
fps = double(fps);
v.FrameRate = fps;
open(v);
%% initialize plots
%nominal trajectory (red)
lv1=1;
Data{lv1} = obj.Data(1,lv1);
Slider{lv1} = patch(Data{lv1}.x1b, Data{lv1}.y1b,'r', 'EdgeAlpha', 1,'FaceAlpha', 1,'EdgeColor', 'r','FaceColor','NONE','LineWidth',0.1);
hold on 
Pusher_a{lv1} = patch(Data{lv1}.X_circle_a,Data{lv1}.Y_circle_a,'r', 'EdgeAlpha', 1,'FaceAlpha', 1, 'EdgeColor', [0,0,1]*0.3,'FaceColor',[1,0,0]*0.5,'LineWidth',0.1);
hold on
Pusher_c{lv1} = patch(Data{lv1}.X_circle_c,Data{lv1}.Y_circle_c,'r', 'EdgeAlpha', 1,'FaceAlpha', 1, 'EdgeColor', [0,0,1]*0.3,'FaceColor',[1,0,0]*0.5,'LineWidth',0.1);
hold on
%actual trajectories (blue)
for lv1=2:obj.NumSim+1
    Data{lv1} = obj.Data(1,lv1); 
    hold on 
    Slider{lv1} = patch(Data{lv1}.x1b, Data{lv1}.y1b,'red', 'EdgeAlpha', 1,'FaceAlpha', 1,'EdgeColor', [0,0,1]*0.3,'FaceColor','NONE','LineWidth',3.0);
    hold on 
    Pusher_a{lv1} = patch(Data{lv1}.X_circle_a,Data{lv1}.Y_circle_a,'r', 'EdgeAlpha', 1,'FaceAlpha', 1, 'EdgeColor', [0,0,1]*0.3,'FaceColor',[1,0,0]*0.5,'LineWidth',0.1);
    hold on
    Pusher_c{lv1} = patch(Data{lv1}.X_circle_c,Data{lv1}.Y_circle_c,'r', 'EdgeAlpha', 1,'FaceAlpha', 1, 'EdgeColor', [0,0,1]*0.3,'FaceColor',[1,0,0]*0.5,'LineWidth',0.1);
    hold on
end
%% update figures
%nominal trajectory (red)
for i1=1:accFactor:length(obj.t)
    lv1=1;
    %update data
    Data{lv1} = obj.Data(i1,lv1);
    %update figures
    Slider{lv1}.XData = Data{lv1}.x1b;
    Slider{lv1}.YData = Data{lv1}.y1b;
    Pusher_a{lv1}.XData = Data{lv1}.X_circle_a;
    Pusher_a{lv1}.YData = Data{lv1}.Y_circle_a;
    Pusher_c{lv1}.XData = Data{lv1}.X_circle_c;
    Pusher_c{lv1}.YData = Data{lv1}.Y_circle_c;
%     %actual trajectories (blue)
    for lv1=2:obj.NumSim+1
        lv1
        Data{lv1} = obj.Data(i1,lv1);
        Slider{lv1}.XData = Data{lv1}.x1b;
        Slider{lv1}.YData = Data{lv1}.y1b;
        Pusher_a{lv1}.XData = Data{lv1}.X_circle_a;
        Pusher_a{lv1}.YData = Data{lv1}.Y_circle_a;
        Pusher_c{lv1}.XData = Data{lv1}.X_circle_c;
        Pusher_c{lv1}.YData = Data{lv1}.Y_circle_c;
        Slider_thin{lv1} = patch(Data{lv1}.x1b, Data{lv1}.y1b,'red', 'FaceAlpha', .2,'EdgeAlpha', .2,'EdgeColor', [0,0,1]*0.3,'FaceColor','NONE','LineWidth',0.1);
        hold on 
        Pusher_thin_a{lv1} = patch(Data{lv1}.X_circle_a,Data{lv1}.Y_circle_a,'red', 'EdgeAlpha', .2,'FaceAlpha', .2, 'EdgeColor', [0,0,1]*0.3,'FaceColor',[1,0,0]*0.5,'LineWidth',0.1);
        hold on 
        Pusher_thin_c{lv1} = patch(Data{lv1}.X_circle_c,Data{lv1}.Y_circle_c,'red', 'EdgeAlpha', .2,'FaceAlpha', .2, 'EdgeColor', [0,0,1]*0.3,'FaceColor',[1,0,0]*0.5,'LineWidth',0.1);
    end
    %update and save frame
    frame = getframe(obj.Ani);
    writeVideo(v,frame);
end
close(v);
end

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
    thetap = x_state(ite, 6);
    u1 = u_state(ite, 1);
    u2 = u_state(ite, 2);
    u3 = u_state(ite, 3);
    %Rotation matrix
    Rb = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    Rp = [cos(thetap) -sin(thetap); sin(thetap) cos(thetap)];
    %kinematics
    Cbi = Helper.C3_2d(theta);
    Cpi = Helper.C3_2d(thetap);
    %position vectors (inertial frame)
    ribi = [x;y];
    ripi = [xp;yp];
    ripb = ripi-ribi;
    %position vectors (other frames)
    rbbi = Cbi*ribi;
    rpap = [0;obj.d];
    rpcp = [0;-obj.d];
    riap = Cpi'*rpap;
    ricp = Cpi'*rpcp;
    rbpb = Cbi*ripb;
    rbpi = Cbi*ripi;
    rbap = Cbi*riap;
    rbcp = Cbi*ricp;
    %find distances
    rbab = rbpi + rbap - rbbi;
    rbcb = rbpi + rbcp - rbbi;
    rbai = rbpi + rbap;
    rbci = rbpi + rbcp;
    rbzi = rbbi-obj.a/2;
    rbaz = rbai-rbzi;
    rbcz = rbci-rbzi;
    npa_b = [1;0];
    npc_b = [1;0];
    rb{1} = rbab;
    rb{2} = rbcb;
    
    %% Calculate distances
    d{1} = -rbaz'*npa_b;
    d{2} = -rbcz'*npc_b;
    %% contact boolean variables
    contact = {};
    for lv1=1:2
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
    %Point a (pusher)
    posa = ripi + Rp*(rpap) + 1*Rb*([-Radius;0]);
    X_circle_a = X+posa(1);
    Y_circle_a = Y+posa(2);
    %Point c (pusher)
    posc = ripi + Rp*(rpcp) + 1*Rb*([-Radius;0]);
    X_circle_c = X+posc(1);
    Y_circle_c = Y+posc(2);
    %Compute velocity vectors
    vbpi = [u1;u2];
    dthetap = u3;
    vipi = Cbi'*vbpi;
    %Set ouputs
    AniOutput.x1b = x1b; 
    AniOutput.y1b = y1b; 
    AniOutput.X_circle_a = X_circle_a;
    AniOutput.Y_circle_a = Y_circle_a;
    AniOutput.X_circle_c = X_circle_c;
    AniOutput.Y_circle_c = Y_circle_c;
    AniOutput.posa = posa;
    AniOutput.posc = posc;
    AniOutput.contact = contact;
end

end

methods
end
    
end

