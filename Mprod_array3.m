% FUNCTION
% Mprod_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of matrix product 
%
% SYNOPSIS
% C = Mprod_array3(A,B)
%
% INPUT
% A, B (i.e., matrices) 
%
% OUTPUT
% C (i.e., matrix)
%
% DESCRIPTION
% Computation, for all frames (i.e., in 3rd dimension, cf. data structure
% in user guide), of the product of two matrices of compatible sizes
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

function C = Mprod_array3(A,B)

% Dimensions
l1 = size(A,1); % Number of lines
c1 = size(A,2); % Number of columns
l2 = size(B,1); % Number of lines
c2 = size(B,2); % Number of columns
n = size(A,3); % Number of frames

% Initilization
C = [];

if l1 == 1 == c1 == 1 % A is scalar
        
    % Element by element product (1*1*n)
        for i = 1:l2
            for j = 1:c2
                C(i,j,:) = A(1,1,:).*B(i,j,:);
            end
        end
        
elseif l2 ~= c1
    
    % Display
    disp('A and B must be of compatible size')
    
else % A and B are matrices

    % transpose = permute ( , [2,1,3])
    At = permute(A,[2,1,3]);
    
    % Element by element product (1*1*n)
    for j =1:c2
        for i = 1:l1
            if l2 == 1 == c1 == 1 % Scalar product
                C(i,j,:) = At(:,i,:).*B(:,j,:);
            else % Dot product (of vectors)
                C(i,j,:) = dot(At(:,i,:),B(:,j,:));
            end
        end
    end
    
end
