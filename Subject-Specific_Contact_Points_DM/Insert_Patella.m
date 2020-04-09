%  FUNCTION
% Insert_Patella.m
%__________________________________________________________________________
%
% PURPOSE
% Insert patella as segment 4
%
% SYNOPSIS
% Segment = Insert_Patella.m(Segment,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model to pass number of frames
%
% OUTPUT
% Segment (cf. data structure in user guide)
%
% DESCRIPTION
% Index incrementation for the thigh and pelvis segments
% Definition of the parameter Q of patella segment
%
% REFERENCES
% N. Sancisi and V. Parenti-Castelli, A new kinematic model of the passive
% motion of the knee inclusive of the patella. Journal of mechanisms and
% robotics 2011, 3(4): 041003.
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Vnorm_array3.m
% Mprod_array3.m
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
% Modify_Segment.m renamed Insert_Patella.m
%__________________________________________________________________________

function Segment = Insert_Patella(Segment,Model)

% Number of frames
n = Model.Informed.n;


%% -------------------------------------------------------------------------
% Model parameters
% -------------------------------------------------------------------------

Segment(4).Informed.rM = []; % No associated marker
Segment(4).Informed.m = 0.0975; % According to Delp's model (no scaling)
Segment(4).Informed.rCs = [0;0;0]; % rCs = rP4 (origin of the SCS)
Segment(4).Informed.Is = zeros(3,3); % Patella is considered as a ponctual mass

invB5 = inv(Segment(5).Informed.B);
nZ5 = invB5(:,3);
% Interpolation matrix
NVZ5 =[nZ5(1,1)*eye(3),...
    (nZ5(2,1))*eye(3), ...
    - nZ5(2,1)*eye(3), ...
    nZ5(3,1)*eye(3)];

% Hinge in patella
rVs14 = Segment(4).Generic.rVs(:,1)*Segment(4).Informed.Scale;

% Patella length
Segment(4).Informed.L = Segment(4).Generic.L*Segment(4).Informed.Scale;

% Hinge in femur
rVs65 = Segment(5).Generic.rVs(:,6)*Segment(5).Informed.Scale;

% Point of the hinge axis in femur
Segment(5).Informed.nV(:,6) =  invB5*rVs65; % Virtual marker 6
% Interpolation matrix
NV65 = [Segment(5).Informed.nV(1,6)*eye(3),...
    (1 + Segment(5).Informed.nV(2,6))*eye(3),...
    - Segment(5).Informed.nV(2,6)*eye(3), ...
    Segment(5).Informed.nV(3,6)*eye(3)];


%% -------------------------------------------------------------------------
% Definition of segment parameters Q and corresponding geometry
% -------------------------------------------------------------------------

u4 = Vnorm_array3((Segment(3).Informed.Q(1:3,1,:) + ...
    Segment(5).Informed.Q(1:3,1,:))/2); % Mean of thigh and shank u axes
v4 = Vnorm_array3(cross(Mprod_array3(repmat(NVZ5,[1,1,n]), ...
    Segment(5).Informed.Q),u4)); % Z axis of femur SCS in the NSCS
w4 = Vnorm_array3(cross(u4,v4));
rP4 = Mprod_array3(repmat(NV65,[1,1,n]),Segment(5).Informed.Q)...
    - rVs14(1,1)*u4 ...
    - rVs14(2,1)*v4; % Distance from the point of the hinge axis to the origin of patella SCS
rD4 = rP4 - Segment(4).Informed.L*v4;
Segment(4).Informed.Q = [u4;rP4;rD4;w4];

% Corresponding geometry
Segment(4).Informed.alpha = Segment(4).Generic.alpha;
Segment(4).Informed.beta = Segment(4).Generic.beta;
Segment(4).Informed.gamma = Segment(4).Generic.gamma;
Segment(4).Informed.B = [1, ...
    Segment(4).Informed.L*cosd(Segment(4).Informed.gamma), ...
    cosd(Segment(4).Informed.beta); ...
    0, ...
    Segment(4).Informed.L*sind(Segment(4).Informed.gamma), ...
    (cosd(Segment(4).Informed.alpha) - ...
    cosd(Segment(4).Informed.beta)*cosd(Segment(4).Informed.gamma))/ ...
    sind(Segment(4).Informed.gamma); ...
    0, ...
    0, ...
    sqrt(1 - cosd(Segment(4).Informed.beta)^2 - ...
    ((cosd(Segment(4).Informed.alpha) - ...
    cosd(Segment(4).Informed.beta)*cosd(Segment(4).Informed.gamma))/ ...
    sind(Segment(4).Informed.gamma))^2)];
Segment(4).Informed.T = Q2Tuv_array3(Segment(4).Informed.Q);


