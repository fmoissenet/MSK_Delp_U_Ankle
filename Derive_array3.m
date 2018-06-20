% FUNCTION
% Derive_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of time derivative of matrix or vector 
%
% SYNOPSIS
% dVdt = Derive_array3(V,dt)
% dMdt = Derive_array3(M,dt)
%
% INPUT
% V (i.e., vector) or M (i.e., matrix) 
% dt (i.e., sampling time that is to say dt = 1/f)
%
% OUTPUT
% dVdt or dMdt (i.e., vector or matrix)
%
% DESCRIPTION
% Gradient approximation in the 3rd dimension of matrix or vector 
% (i.e., all frames, cf. data structure in user guide)
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

function dVdt = Derive_array3(V,dt)

% Trick for passing both versions R14 and R2007
if size(V,2) == 1
    [gr1,gr2,dVdt] = gradient(repmat(V,[1,size(V,1),1]),dt);
    dVdt = dVdt(1:size(V,1),1,:);
else
    % 3rd direction gradient
    % with time sampling (dt)
    [gr1,gr2,dVdt] = gradient(V,dt);
end
