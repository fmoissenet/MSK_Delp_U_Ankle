% FUNCTION
% Vnop_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Non-orthogonal projection on three basis vectors
%
% SYNOPSIS
% Vnop = Vnop_array3(V,e1,e2,e3)
%
% INPUT
% V (i.e., vector)
% e1, e2, e3 (i.e., basis vectors)
%
% OUTPUT
% Vnop (i.e., vector components)
%
% DESCRIPTION
% Computation, for all frames (i.e., in 3rd dimension, cf. data structure
% in user guide), of the vector projections on three non-orthogonal basis
% vectors
%
% REFERENCES
% L Cheze. Comparison of different calculations of three-dimensional joint
% kinematics from video-based system data. Journal of Biomechanics 2000;
% 33(12):1695-9.
% G Desroches, L Cheze, R Dumas. Expression of joint moment in the joint 
% coordinate system. Journal of Biomechanical Engineering 2010;132(11):
% 114503
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX) 
% None
%
% MATLAB VERSION
% Matlab R2007b
%__________________________________________________________________________
%
% CHANGELOG
% Created by Rapha�l Dumas
% March 2010
%
% Modified by Rapha�l Dumas
% September 2012
% Test of the squared norm for gimbal lock occurence
%__________________________________________________________________________
%
% Copyright (C) 2018  Raphael Dumas, Florent Moissenet
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%__________________________________________________________________________

function Vnop = Vnop_array3(V,e1,e2,e3)

% Projection in a non-orthonormal basis (e1,e2,e3)
% using mixted products
% V = ((e2,e3,V)/(e1,e2,e3))*e1
%      + ((e3,e1,V)/(e1,e2,e3))*e2
%      + ((e1,e2,V)/(e1,e2,e3))*e3
Vnop(1,1,:) = dot(cross(e2,e3),V)./ ...
    dot(cross(e1,e2),e3);
Vnop(2,1,:) = dot(cross(e3,e1),V)./ ...
    dot(cross(e1,e2),e3);
Vnop(3,1,:) = dot(cross(e1,e2),V)./ ...
    dot(cross(e1,e2),e3);

% % Test of the squared norm for gimbal lock occurence
% [ind,~] = find(permute(Vnop(1,1,:).^2 + Vnop(2,1,:).^2 + Vnop(3,1,:).^2 + ...
%     2*(Vnop(1,1,:).*Vnop(2,1,:).*dot(e1,e2) +  ...
%     Vnop(1,1,:).*Vnop(3,1,:).*dot(e1,e3) + ...
%     Vnop(2,1,:).*Vnop(3,1,:).*dot(e2,e3)) - ...
%     sum(V.^2),[3,2,1]) > 1e-3);
% if ~isempty(ind)
%     display(['Warning: gimbal lock occurence at frame(s) ',num2str(ind')])
% end

