% FUNCTION
% Vskew_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of skew matrix of vector
%
% SYNOPSIS
% M = Vskew_array3(V)
%
% INPUT
% V (i.e., vector)
%
% OUTPUT
% M (i.e., matrix)
%
% DESCRIPTION
% Computation, for all frames (i.e., in 3rd dimension, cf. data structure
% in user guide), of the skew matrix from a vector
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

function M =  Vskew_array3(V)

% Number of frame
n = size(V,3); 

% Initialization
M = zeros(3,3,n); % 3*3*n zeros matrix

% Non-zero terms
M(1,2,:) = -V(3,1,:);
M(1,3,:) = V(2,1,:);
M(2,1,:) = V(3,1,:);
M(2,3,:) = -V(1,1,:);
M(3,1,:) = -V(2,1,:);
M(3,2,:) = V(1,1,:);
