% FUNCTION
% Vnorm_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of a unitary vector
%
% SYNOPSIS
% Vn = Vnorm_array3(V)
%
% INPUT
% V (i.e., vector)
%
% OUTPUT
% Vn (i.e., vector)
%
% DESCRIPTION
% Computation, for all frames (i.e., in 3rd dimension, cf. data structure
% in user guide), of the unitary vector from a vector
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

function Vn = Vnorm_array3(V)

% Norm in dimension 1*n
N = sqrt(sum(V.^2));

% Number of line
L = size(V,1);
% Matrix for normalization in dimension (L*L*n)
for i = 1:L
    In(i,i,:) = 1/N;
end
Vn = Mprod_array3(In,V);
