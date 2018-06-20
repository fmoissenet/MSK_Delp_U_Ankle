% FUNCTION
% Q2Twu_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Alternative computation of homogenous matrix from generalized coordinates
%
% SYNOPSIS
% T = Q2Twu_array3(Q)
%
% INPUT
% Q (i.e., generalized coordinates)
%
% OUTPUT
% T (i.e., homogenous matrix)
%
% DESCRIPTION
% Computation, for all frames (i.e., in 3rd dimension, cf. data structure
% in user guide), of the homogenous matrix (T) from generalized coordinates
% (Q) with axis correspondence w = Z and origin at endpoint D
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX) 
% Mprod_array3.m
% Minv_array3.m
%
% MATLAB VERSION
% Matlab R2007b
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas
% March 2010
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

function T = Q2Twu_array3(Q)

L = sqrt(sum((Q(4:6,1,:) - Q(7:9,1,:)).^2));
% Alpha angle between (rP - rD) and w
a = acos(dot(Q(4:6,1,:) - Q(7:9,1,:),Q(10:12,1,:))./...
    sqrt(sum((Q(4:6,1,:) - Q(7:9,1,:)).^2)))*180/pi;
% Beta angle between u and w
b = acos(dot(Q(1:3,1,:),Q(10:12,1,:)))*180/pi;
% Gamma angle between u and (rP - rD)
c = acos(dot(Q(1:3,1,:),Q(4:6,1,:) - Q(7:9,1,:))./...
    sqrt(sum((Q(4:6,1,:) - Q(7:9,1,:)).^2)))*180/pi;

% B matrix from (u rP-rD w) to (X Y Z)
B(1,1,:) = sind(b);
B(1,2,:) = (L.*cosd(c)-L.*cosd(a).*cosd(b))./sind(b);
B(2,2,:) = L.*sqrt(1-(cosd(a).^2)-((cosd(c)-cosd(a).*cosd(b))./sind(b)).^2);
B(3,1,:) = cosd(b);
B(3,2,:) = L.*cosd(a);
B(3,3,:) = 1;

% [X Y Z rP] = [[u rP-rD w]*B-1 rD]
T = [Mprod_array3([Q(1:3,1,:),Q(4:6,1,:) - Q(7:9,1,:),Q(10:12,1,:)],...
    Minv_array3(B)),Q(7:9,1,:)];
T(4,4,:) = 1;

% % Alternatively
% % Orthonormal axis
% Z = Q(10:12,1,:); % Z = w
% Y = Vnorm_array3(cross(Z,Q(1:3,1,:))); % Z x u
% X = Vnorm_array3(cross(Y,Z));
% 
% % Homogenous matrix
% T = [X, Y, Z, Q(7:9,1,:)]; 
% T(4,4,:) = 1;
