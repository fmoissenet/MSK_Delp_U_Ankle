% FUNCTION
% Minv_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of matrix inverse 
%
% SYNOPSIS
% B = Minv_array3(A)
%
% INPUT
% A (i.e., matrix) 
%
% OUTPUT
% B (i.e., matrix)
%
% DESCRIPTION
% Computation of the inverse of a square matrix made for all frames (i.e.,
% in 3rd dimension, cf. data structure in user guide) in the case of 
% rotation matrix (R)
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

function B = Minv_array3(A)

if size(A,2) == 3 & ...
        mean(dot(cross(A(:,1,:),A(:,2,:)),A(:,3,:)),3) == 1
    % Determinant 3x3 is mixt product

    % Transposition of rotation matrix
    % with transpose = permute( ,[2,1,3])
    B = permute(A,[2,1,3]);

else % Clasical inversion frame by frame
    for n = 1:size(A,3)
        B(:,:,n) = inv(A(:,:,n));
    end
end


