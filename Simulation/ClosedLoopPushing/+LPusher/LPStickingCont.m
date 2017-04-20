classdef LPStickingCont < LPusher.LPController
    %LPSTICKING Summary of this class goes here
    %   Detailed explanation goes here
    
properties
    name = 'LPSticking';
end

methods
        
function [LocalAin, Localbin, LocalAeq, Localbeq] = buildModeDepConstraints(obj, Opt, index)
    % Sticking Constraint
    numConstraints = 1;
    uIndex = 5;
    Bleft  = 1;
    cLeft  = 0;
    Bright = 0;
    cRight = 0;
    [LocalAin, Localbin, LocalAeq, Localbeq] = obj.addPusherConstraints(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '==', 'star', Opt);
    clear Bleft cLeft Bright cRight
end
end
    
end