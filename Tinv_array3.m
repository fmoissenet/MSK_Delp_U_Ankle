% FUNCTION
% Tinv_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of homogenous matrix inverse 
%
% SYNOPSIS
% B = Tinv_array3(A)
%
% INPUT
% A (i.e., homogenous matrix) 
%
% OUTPUT
% B (i.e., homogenous matrix)
%
% DESCRIPTION
% Computation, for all frames (i.e., in 3rd dimension, cf. data structure
% in user guide), of the inverse of an homogenous matrix (T)
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX) 
% Mprod_array3.m
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

function B = Tinv_array3(A)

% Initialization
A1 = [];
A2 = [];

% First matrix
% Transposition of rotational part
% with transpose = permute( ,[2,1,3])
A1(1:3,1:3,:) = permute(A(1:3,1:3,:),[2,1,3]);
A1(4,4,:) = 1;

% Second matrix
% Opposite sign of the translation part
A2(1:4,4,:) = - A(1:4,4,:);
A2(1,1,:) = 1;
A2(2,2,:) = 1;
A2(3,3,:) = 1;
A2(4,4,:) = 1;

% Matrix product
B = Mprod_array3(A1,A2);
