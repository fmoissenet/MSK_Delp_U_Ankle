% FUNCTION
% Compute_J.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of the pseudo-inertia matrix from the segment mass, 
% inertia and CoM position
%
% SYNOPSIS
% Segment = Compute_J(Segment,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model to pass number of frames
%
% OUTPUT
% Segment (cf. data structure in user guide)
%
% DESCRIPTION
% Computation of the inertia matrix at proximal endpoint in (u,rP-rD,w)
% using generalised paralell axis theorem 
%
% REFERENCES
% 
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Mprod_array3.m
% Minv_array3.m
%
% MATLAB VERSION
% Matlab R2020a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas, Florent Moissent, Edouard Jouan
% September 2012
%
% Modified by Raphael Dumas
% March 2020
% Generic vs. Informed structures (Segment, Joint, Model)
% n,f,fc as Model.Informed fields
%__________________________________________________________________________


function Segment = Compute_J(Segment,Model)

% Number of frames
n = Model.Informed.n;

% Initialisation
E33 = eye(3,3);

for i = 2:5 % From Foot (i = 2) to Thigh (i = 5)
    
        % J = invB*(Is + m*((rCs'*rCs)*E33 - rCs*rCs')*invB'
        Segment(i).Informed.J = ...
            Mprod_array3(Mprod_array3(Minv_array3(repmat(Segment(i).Informed.B,[1,1,n])),...
            repmat(Segment(i).Informed.Is + ...
            Segment(i).Informed.m*((Segment(i).Informed.rCs'*Segment(i).Informed.rCs)*E33 - ...
            Segment(i).Informed.rCs*Segment(i).Informed.rCs'),[1,1,n])),...
            permute(Minv_array3(repmat(Segment(i).Informed.B,[1,1,n])),[2,1,3]));
        % with transpose = permute( ,[2,1,3])
  
end