% FUNCTION
% R2mobileZYX_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of Euler angles from rotation matrix (with ZYX mobile sequence
% for joint kinematics) 
%
% SYNOPSIS
% Joint_Euler_Angles = R2mobileZYX_array3(R)
%
% INPUT
% R (i.e., rotation matrix) 
%
% OUTPUT
% Joint_Euler_Angles (i.e., tetha1, tetha2, tetha3, in line)
%
% DESCRIPTION
% Computation, for all frames (i.e., in 3rd dimension, cf. data structure
% in user guide), of the Euler angles (tetha1, tetha2, tetha3) from the
% rotation matrix (R) using a sequence of mobile axes ZYX
%
% REFERENCE
% L Cheze, R Dumas, JJ Comtet, C Rumelhart, M Fayet. A joint coordinate 
% system proposal for the study of the trapeziometacarpal joint kinematics
% Computer Methods in Biomechanics and Biomedical Engineering 2009;12(3): 
% 277-82
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

function Joint_Euler_Angles = R2mobileZYX_array3(R)

% Tetha1: Flexion-Extension (about Z proximal SCS axis)
Joint_Euler_Angles(1,1,:) = atan2(R(2,1,:),R(1,1,:));
 % Tetha2: Internal-External Rotation (about Y floating axis)
Joint_Euler_Angles(1,2,:) = asin(-R(3,1,:));
% Tetha3: Abduction-Adduction (about X distal SCS axis) 
Joint_Euler_Angles(1,3,:) = atan2(R(3,2,:),R(3,3,:)); 
