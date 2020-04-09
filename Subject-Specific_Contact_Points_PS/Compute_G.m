% FUNCTION
% Compute_G.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of the generalised mass matrix
%
% SYNOPSIS
% [Segment,Model] = Compute_G(Segment,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% OUTPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% DESCRIPTION
% Computation of the segment and model mass matrices 
%
% REFERENCES
% R Dumas, L Cheze. 3D inverse dynamics in non-orthonormal segment 
% coordinate system. Medical & Biological Engineering & Computing 2007;
% 45(3):315-22
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX)
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


function [Segment,Model] = Compute_G(Segment,Model)

% Number of frames
n = Model.Informed.n;

% Initialmisation
E33 = eye(3,3);
G = zeros(48,48);

for i = 2:5 % From Foot (i = 2) to Thigh (i = 5)
    
    % CoM Coordinates in generalized coordinates (u,rP-rD,w)
    % ---------------------------------------------------------------------
    Segment(i).Informed.nC = inv(Segment(i).Informed.B)*Segment(i).Informed.rCs;
    
    % 1st line
    % ---------------------------------------------------------------------
    Segment(i).Informed.G(1:3,1:3) = Segment(i).Informed.J(1,1)*E33;
    Segment(i).Informed.G(1:3,4:6) = (Segment(i).Informed.m*Segment(i).Informed.nC(1,1) + ...
        Segment(i).Informed.J(1,2))*E33;
    Segment(i).Informed.G(1:3,7:9) = - Segment(i).Informed.J(1,2)*E33;
    Segment(i).Informed.G(1:3,10:12) = Segment(i).Informed.J(1,3)*E33;
    
    % 2nd line
    % ---------------------------------------------------------------------
    Segment(i).Informed.G(4:6,1:3) = (Segment(i).Informed.m*Segment(i).Informed.nC(1,1) + ...
        Segment(i).Informed.J(1,2))*E33;
    Segment(i).Informed.G(4:6,4:6) = (Segment(i).Informed.m + ...
        2*Segment(i).Informed.m*Segment(i).Informed.nC(2,1) + Segment(i).Informed.J(2,2))*E33;
    Segment(i).Informed.G(4:6,7:9) = (- Segment(i).Informed.m*Segment(i).Informed.nC(2,1) - ...
        Segment(i).Informed.J(2,2))*E33;
    Segment(i).Informed.G(4:6,10:12) = (Segment(i).Informed.m*Segment(i).Informed.nC(3,1) + ...
        Segment(i).Informed.J(2,3))*E33;
    
    % 3rd line
    % ---------------------------------------------------------------------
    Segment(i).Informed.G(7:9,1:3) = - Segment(i).Informed.J(1,2)*E33;
    Segment(i).Informed.G(7:9,4:6) = (- Segment(i).Informed.m*Segment(i).Informed.nC(2,1) - ...
        Segment(i).Informed.J(2,2))*E33;
    Segment(i).Informed.G(7:9,7:9) = Segment(i).Informed.J(2,2)*E33;
    Segment(i).Informed.G(7:9,10:12) = - Segment(i).Informed.J(2,3)*E33;
    
    % 4th line
    % ---------------------------------------------------------------------
    Segment(i).Informed.G(10:12,1:3) = Segment(i).Informed.J(1,3)*E33;
    Segment(i).Informed.G(10:12,4:6) = (Segment(i).Informed.m*Segment(i).Informed.nC(3,1) + ...
        Segment(i).Informed.J(2,3))*E33;
    Segment(i).Informed.G(10:12,7:9) = - Segment(i).Informed.J(2,3)*E33;
    Segment(i).Informed.G(10:12,10:12) = Segment(i).Informed.J(3,3)*E33;
    
    % Creation of the generalized mass matrix (48x48)
    % ---------------------------------------------------------------------
    s = 12*(i-2);
    for j = 1:12
        for k = 1:12
            G(s+j,s+k) = Segment(i).Informed.G(j,k);
        end
    end
    
end

% Generalized mass matrix 
Model.Informed.G = repmat(G,[1,1,n]);
