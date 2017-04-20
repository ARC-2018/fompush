classdef LPUpCont < LPusher.LPController
    %LPUP Summary of this class goes here
    %   Detailed explanation goes here
    
properties
    name = 'LPDown';
end

methods
function [Ai, bi, Ae, be] = buildModeDepConstraints(obj, Opt, index)
    %Sliding left constraint
    numConstraints = 1;
    uIndex = [5];
    Bleft  = [1];
    cLeft  = [0];
    Bright = [0];
    cRight = [0];
    [Ai, bi, Ae, be] = obj.addPusherConstraints2(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '>', 'star', Opt);
    clear Bleft cLeft Bright cRight
    %friction cone edge constraints
    numConstraints = 1;
    uIndex = [1,3];
    Bleft  = [0 1];
    cLeft  = [0];
    Bright = [obj.nu_p 0];
    cRight = [0];
    [A1,b1,A2,b2] = obj.addPusherConstraints2(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '==', 'star', Opt);
    Ai = [Ai; A1];
    bi = [bi; b1];
    Ae = [Ae; A2];
    be = [be; b2];
    clear Bleft cLeft Bright cRight
    %friction cone edge constraints
    numConstraints = 1;
    uIndex = [2,4];
    Bleft  = [0 1];
    cLeft  = [0];
    Bright = [obj.nu_p 0];
    cRight = [0];
    [A1,b1,A2,b2] = obj.addPusherConstraints2(Bleft, Bright, cLeft, cRight, numConstraints, uIndex, index, '==', 'star', Opt);
    Ai = [Ai; A1];
    bi = [bi; b1];
    Ae = [Ae; A2];
    be = [be; b2];
    clear Bleft cLeft Bright cRight
end
        
end
    
end

