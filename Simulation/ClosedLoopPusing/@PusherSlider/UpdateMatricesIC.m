function obj = UpdateMatricesIC(obj,counter,delta_x,Family)
    rx = -obj.a/2;
    ry = delta_x(4);
    gamma_top     = ( obj.nu_p*obj.c^2 - rx*ry + obj.nu_p*rx^2)/(obj.c^2 + ry^2 - obj.nu_p*rx*ry);
    gamma_bottom  = (-obj.nu_p*obj.c^2 - rx*ry - obj.nu_p*rx^2)/(obj.c^2 + ry^2 + obj.nu_p*rx*ry);
    %% Add nonlinear dynamic constraint
    B = obj.B_Nonlinear{Family}(ry,delta_x(3));
    F = B*obj.u_star - obj.f_star{Family}(0,0,obj.u_star(1),obj.u_star(2));
    F_tilde = delta_x + obj.h_opt*F;
    B_tilde = obj.h_opt*B;
    %Add constraint (Modify existing dummy constraint)
    Aeq = zeros(obj.num_vars, obj.Opt.FOM{Family}.nv);
    Aeq(:,obj.Opt.FOM{Family}.vars.x.i(1:obj.num_vars,1))   = eye(obj.num_vars);
    Aeq(:,obj.Opt.FOM{Family}.vars.u.i(1:obj.num_inputs,1)) = -B_tilde;
    beq = F_tilde;
    Ain1 = Aeq;
    bin1 = beq;
    Ain2 = -Aeq;
    bin2 = -beq;
    Ain = [Ain1;Ain2];
    bin = [bin1;bin2];
    obj.Opt.FOM{Family} = obj.Opt.FOM{Family}.addLinearConstraints(Ain, bin, [], []);
%     obj.Opt.FOM{Family} = obj.Opt.FOM{Family}.addLinearConstraints(Ain2, bin2, [], []);
    clear Aeq beq
    %% Add nonlinear motion cone constraint
    epsilon = 0.005;
    NconstMC = 1;
    if counter==1
            %Sticking MC (2 constraints)
            D = [-gamma_top 1];
            bin = -obj.u_star(2) + gamma_top*obj.u_star(1);
            Ain = zeros(NconstMC, obj.Opt.FOM{Family}.nv);
            Ain(:,obj.Opt.FOM{Family}.vars.u.i(1:2,1)) = D;
%             obj.Opt.FOM{Family}.A(9,:) = Ain;
%             obj.Opt.FOM{Family}.b(9,1) =  bin;
            obj.Opt.FOM{Family} = obj.Opt.FOM{Family}.addLinearConstraints(Ain, bin, [], []);   
            clear D E Ain bin
            D = [gamma_bottom -1];
            bin = obj.u_star(2) - gamma_bottom*obj.u_star(1);
            Ain = zeros(NconstMC, obj.Opt.FOM{Family}.nv);
            Ain(:,obj.Opt.FOM{Family}.vars.u.i(1:2,1)) = D;
%             obj.Opt.FOM{Family}.A(10,:) = Ain;
%             obj.Opt.FOM{Family}.b(10,1) =  bin;
            obj.Opt.FOM{Family} = obj.Opt.FOM{Family}.addLinearConstraints(Ain, bin, [], []); 
            clear D E Ain bin
        elseif counter==2
            %Sliding up MC (1 constraint)
            D = [gamma_top -1];
            bin = obj.u_star(2) - gamma_top*obj.u_star(1)-epsilon;
            Ain = zeros(NconstMC, obj.Opt.FOM{Family}.nv);
            Ain(:,obj.Opt.FOM{Family}.vars.u.i(1:2,1)) = D;
%             obj.Opt.FOM{Family}.A(9,:) = Ain;
%             obj.Opt.FOM{Family}.b(9,1) =  bin;
%             obj.Opt.FOM{Family}.A(10,:) = Ain*0;
%             obj.Opt.FOM{Family}.b(10,1) =  0;
            obj.Opt.FOM{Family} = obj.Opt.FOM{Family}.addLinearConstraints(Ain, bin, [], []);   
            clear D E Ain bin
        else
            %Sliding down MC (1 constraint)
            D = [-gamma_bottom 1];
            bin = -obj.u_star(2) + gamma_bottom*obj.u_star(1)-epsilon;
            Ain = zeros(NconstMC, obj.Opt.FOM{Family}.nv);
            Ain(:,obj.Opt.FOM{Family}.vars.u.i(1:2,1)) = D;
%             obj.Opt.FOM{Family}.A(9,:) = Ain;
%             obj.Opt.FOM{Family}.b(9,1) =  bin;
%             obj.Opt.FOM{Family}.A(10,:) = Ain*0;
%             obj.Opt.FOM{Family}.b(10,1) =  0;
            obj.Opt.FOM{Family} = obj.Opt.FOM{Family}.addLinearConstraints(Ain, bin, [], []);   
            clear D E Ain bin
    end
%     obj.Opt.FOM{Family}.A
%     obj.Opt.FOM{Family}.b
end