classdef LPStickingCont < LPusher.LPController
    %LPSTICKING Summary of this class goes here
    %   Detailed explanation goes here
    
properties
    name = 'LPSticking';
end

methods
        
function [Ai, bi, Ae, be] = buildModeDepConstraints(obj, Opt, index)
    % Sticking Constraint
    numConstraints = 1;
    uIndex = 5;
    Bleft  = 1;
    cLeft  = 0;
    Bright = 0;
    cRight = 0;
    [Ai, bi, Ae, be] = obj.addPusherConstraints2(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '==', 'star', Opt);
    clear Bleft cLeft Bright cRight
end
end
    
end