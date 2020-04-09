% FUNCTION
% Compute_P.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of the generalised weight
%
% SYNOPSIS
% Model = Compute_P(Segment,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% OUTPUT
% Model (cf. data structure in user guide)
%
% DESCRIPTION
%
% REFERENCES
% R Dumas, L Cheze. 3D inverse dynamics in non-orthonormal segment 
% coordinate system. Medical & Biological Engineering & Computing 2007;
% 45(3):315-22
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
%
% MATLAB VERSION
% Matlab R2020a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas, Florent Moissent, Edouard Jouan
% September 2012
% Modified by Raphael Dumas
%
% March 2020
% Generic vs. Informed structures (Segment, Joint, Model)
% n,f,fc as Model.Informed fields
%__________________________________________________________________________


function Model = Compute_P(Segment,Model)

% Number of frames
n = Model.Informed.n;

% Initialisation
E33 = eye(3,3);
P = [];

for i = 2:5 % From Foot (i = 2) to Thigh (i = 5) 
P = [P; ...
    [Segment(i).Informed.nC(1,1)*E33,...
    (1 + Segment(i).Informed.nC(2,1))*E33,...
    - Segment(i).Informed.nC(2,1)*E33,...
    Segment(i).Informed.nC(3,1)*E33]'*Segment(i).Informed.m*[0;-9.81;0]];
end

% Generalised weight
Model.Informed.P = repmat(P,[1,1,n]);
