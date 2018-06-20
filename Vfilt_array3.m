% FUNCTION
% Vfilt_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Filtering of vector 
%
% SYNOPSIS
% Vf = Vfilt_array3(V,f,fc)
%
% INPUT
% V (i.e., vector) 
% f (i.e., sampling frequency)
% fc (i.e., cut frequency)
%
% OUTPUT
% Vf (i.e., vector)
%
% DESCRIPTION
% Filtering, along with the 3rd dimension (i.e., all frames, cf. data
% structure in user guide), of the vector components by a 4th order
% Butterworth with special attention when the vector is a column of an
% homogenous matrix
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX) 
% None
% 
% MATLAB VERSION
% Matlab R2013a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas
% March 2010
%
% Modified by Raphael Dumas
% September 2012
% Filtering line by line (any number of line)
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

function Vf = Vfilt_array3(V,f,fc)

% Butterworth
[af,bf] = butter(4,fc./(f/2));

% Initialisation
Vf = [];
% Filtering line by line
for c = 1:size(V,1)
    Vc = filtfilt(af,bf,permute(V(c,1,:),[3,1,2]));
    Vf = [Vf;permute(Vc,[3,2,1])];
end

