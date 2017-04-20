classdef LPDownCont < LPusher.LPController
    %LPDOWN Summary of this class goes here
    %   Detailed explanation goes here
    
properties
    name = 'LPDown';
end

methods
function [Opt] = buildModeDepConstraints(obj, Opt, index)
%   Sliding left constraint
    numConstraints = 1;
    uIndex = [5];
    Bleft  = [1];
    cLeft  = [0];
    Bright = [0];
    cRight = [0];
    [Opt] = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '<=', 'star', Opt);
    clear Bleft cLeft Bright cRight
    %friction cone edge constraints
    numConstraints = 1;
    uIndex = [1,3];
    Bleft  = [0 1];
    cLeft  = [0];
    Bright = [-obj.nu_p 0];
    cRight = [0];
    [Opt] = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '==', 'star', Opt);
    clear Bleft cLeft Bright cRight
%   friction cone edge constraints
    numConstraints = 1;
    uIndex = [2,4];
    Bleft  = [0 1];
    cLeft  = [0];
    Bright = [-obj.nu_p 0];
    cRight = [0];
    [Opt] = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '==', 'star', Opt);
    clear Bleft cLeft Bright cRight
%     LocalAin = [Ain1; Ain2; Ain3];
%     Localbin = [bin1; bin2; bin3];
%     LocalAeq = [Aeq1; Aeq2; Aeq3];
%     Localbeq = [beq1; beq2; beq3];
    
end
        
end

end

