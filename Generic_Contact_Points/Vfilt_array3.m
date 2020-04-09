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

