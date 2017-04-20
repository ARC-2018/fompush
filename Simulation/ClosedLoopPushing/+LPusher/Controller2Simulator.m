function [xs, us] = Controller2Simulator(xc, uc)
%CONTROLLER2SIMULATOR Summary of this function goes here
%   Variables
    friction = LPusher.Friction;
    ribi = xc(1:2);
    theta = xc(3);
    ry = xc(4);
    %Direction cosine matrices
    Cbi = Helper.C3_2d(theta);
    %Convert state to body frame
    rbbi = Cbi * ribi;
    rbpb = [-LPusher.LPusherSlider.a / 2;ry];
    ripb = Cbi'*rbpb;
    ripi = ribi+ripb;
    %Forces
    xs = [ribi;theta;ripi;theta];
    if nargin > 1
        fn = uc(1:2);
        ft = uc(3:4);
        dry = uc(5);
        drbpb = [0;dry];
        %Find rx, ry: In body frame   
        rbap = [0;LPusher.LPusherSlider.d];
        rbcp = [0;-LPusher.LPusherSlider.d];    
        rbab = rbpb + rbap;
        rbcb = rbpb + rbcp;
        %kinematics
        npa_b = [1;0];
        tpa_b = [0;1];
        Tpa_b = [tpa_b];
        Jpa = [1 0 -rbab(2);...
               0 1 rbab(1)];
        %kinematics
        npc_b = [1;0];
        tpc_b = [0;1];
        Tpc_b = [tpc_b];
        Jpc = [1 0 -rbcb(2);...
               0 1 rbcb(1)];
        %build useful matrices pt a
        N{1} = npa_b'*Jpa;
        T{1} = Tpa_b'*Jpa;
        %build useful matrices pt a
        N{2} = npc_b'*Jpc;
        T{2} = Tpc_b'*Jpc;
        %concatenate matrices
        N_tot= [];
        T_tot=[];
        for lv1=1:2
            N_tot = [N_tot;N{lv1}];
            T_tot = [T_tot;T{lv1}];
        end
        %Compute twist (in body frame)
        twist_body_b = friction.A_ls*(N_tot'*fn+T_tot'*ft);
        %convert to inertial frame
        S_ib = [Cbi' [0;0];0 0 1];
        twist_body_i = S_ib*twist_body_b;
        %Compute twist of pusher
        twist_pusher_i = zeros(3,1);
        twist_pusher_i(1:2) = twist_body_i(1:2) + Cbi'*drbpb + Helper.cross3d(twist_body_i(3), ripb);
        twist_pusher_i(3) = twist_body_i(3);
        %output
        us = twist_pusher_i;
    else
        us = nan;
    end
end