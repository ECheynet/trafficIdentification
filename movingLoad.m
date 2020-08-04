function [F] = movingLoad(t,tStart,L,V,x,mass,direction,varargin)
%  [F] = movingLoad(t,tStart,L,V,x,mass,direction) computes the force from
%  moving mass on a beam, discretized with a constant spacing.
% 
% Input:
%  - t: [1 x N] double: time
%  - tStart: [1 x 1] double: Arrival time of the mass
%  - L: [1 x 1] double: length of the beam
%  - x: [1 x M] double: normalized coordinates of the nodes discretizing the
%       beam ( 0 <= x <=1 )
%  - V: [1 x 1] double: speed of the moving mass
%  - mass: [1 x 1] double: mass of the moving mass
%  - Direction: 1 or -1 (either left to right or right to left)
% 
% Output
%  - F: [1 x M]: Force on each nodes of the beam
% 
% Author: E. Cheynet - UiB - 03-08-2020

%%
%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('g',9.81);
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
g = p.Results.g; % acceleration of gravity for earth


%%
F = zeros(1,numel(x));
% time  at which the mass arrives on the beam
if t>=tStart
    pos = V*(t-tStart);
    if direction~=1,            pos = L-pos;end
    
    if pos <= L && pos >0
        [~,indPos]=min(abs(x.*L-pos));
        F(indPos) = -g*mass./median(diff(L.*x));
    end
end

end

