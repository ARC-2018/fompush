classdef Simulation < dynamicprops & PusherSliderSystem
    properties (Constant)

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
        function obj = Simulation(SimName)  
            %% Save data in new folder
            obj.SimName = SimName;
            if ismac()
                tempName = strcat('/Users/Francois/Dropbox (MIT)/Data/',obj.SimName);
            elseif isunix() && ~ismac()
                homeName = getenv('HOME');
                tempName = strcat(homeName, '/Data/cpush/',obj.SimName);
            end
            mkdir(tempName);
            obj.FilePath = tempName;
            obj.FileName = strcat(obj.FilePath,'/',obj.SimName);
            %Nominal Trajectory  #Frank hack
            %Initialize variables
            c = Controller();
            s = Simulator();
            obj.xs_state = zeros(obj.N, s.num_xsStates);
            obj.us_state = zeros(obj.N, s.num_usStates);
            obj.xc_state = zeros(obj.N, c.num_xcStates);
            obj.uc_state = zeros(obj.N, c.num_ucStates);
            obj.xs_star_state= zeros(obj.N, s.num_xsStates);
            obj.us_star_state= zeros(obj.N, s.num_usStates);
            obj.xc_star_state= zeros(obj.N, c.num_xcStates);
            obj.uc_star_state= zeros(obj.N, c.num_ucStates);

        end
    end
end