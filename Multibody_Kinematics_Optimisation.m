% FUNCTION
% Multibody_Kinematics_Optimisation.m
%__________________________________________________________________________
%
% PURPOSE
% Multibody optimisation of the postions of foot, shank, patella, thigh
% and computation of their consistent velocities and accelerations with
% universal ankle joint, parallel tibio-femoral mechanism,
% paralell patello-femoral mechanism and spherial hip joint
%
% SYNOPSIS
% [Segment,Joint] = Multibody_Kinematics_Optimisation(Segment,Joint,f)
%
% INPUT
% Segment (cf. data structure in user guide)
% Joint (cf. data structure in user guide)
% f (frequency)
%
% OUTPUT
% Segment (cf. data structure in user guide)
% Joint (cf. data structure in user guide)

%
% DESCRIPTION
% Computation of Q by minimisation under constraints (by Gauss-Newton)
% Computation of the velocities and accelerations 
%
% Joint models are universal ankle joint, parallel tibio-femoral
% mechanism, paralell patello-femoral mechanism, spherial hip joint
%
% REFERENCES
% F. Moissenet, L. Cheze, R. Dumas, Anatomical kinematic constraints: 
% consequences on musculo-tendon forces and joint reactions. Multibody 
% System Dynamics 2012, 28:125–141.
% S Duprey, L. Cheze, R. Dumas. Influence of joint constraints on lower 
% limb kinematics estimation from skin markers using global optimization.
% Journal of Biomechancis 2010, 19;43(14):2858-62.
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Vnop_array3.m
% Mprod_array3.m
% Derive_array3.m
% Vfilt_array3.m
%
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas, Florent Moissent, Edouard Jouan
% September 2012
%
% Modified by Raphael Dumas
% February 2015
% Jacobian of the kinematic constraints for spherical joints in SCS
%
% Modified by Raphael Dumas
% April 2015
% Change Segment().a into Segment().alpha
% Change Segment().b into Segment().beta
% Change Segment().c into Segment().gamma
%
% Modified by Raphael Dumas
% November 2017
% Universal joint at the ankle according to Delp
%
% Modified by Raphael Dumas
% January 2018
% Renamed to Multibody Kinematics Optimisation
% Position of talus joint centre
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


function [Segment,Joint] = Multibody_Kinematics_Optimisation(Segment,Joint,f)

% Number of frames
n = size(Segment(2).rM,3);
% Initialisation
Joint(1).Kk = [];
fc = 5; % Cut frequency for filtering


%% -------------------------------------------------------------------------
% Model parameters
% -------------------------------------------------------------------------
% Initialisation
Segment(1).L = NaN; % No value for segment 1 (Forceplate)
Segment(1).alpha = NaN; % No value for segment 1 (Forceplate)
Segment(1).beta = NaN; % No value for segment 1 (Forceplate)
Segment(1).gamma = NaN; % No value for segment 1 (Forceplate)
Joint(1).d(1,1:6) = NaN(1,6); % No value for joint 1 (CoP)
Joint(2).d(1,1:6) = NaN(1,6); % No value for joint 2

% Mean segment geometry and markers coordinates
for i = 2:6 % From i = 2 (Foot) to i = 6 (Pelvis)
    
    % Segment length
    Segment(i).L = mean(sqrt(sum((Segment(i).Q(4:6,1,:) - ...
        Segment(i).Q(7:9,1,:)).^2)),3);
    
    % Alpha angle between (rP - rD) and w
    Segment(i).alpha = mean(acosd(dot(Segment(i).Q(4:6,1,:) - ...
        Segment(i).Q(7:9,1,:), Segment(i).Q(10:12,1,:))./...
        sqrt(sum((Segment(i).Q(4:6,1,:) - ...
        Segment(i).Q(7:9,1,:)).^2))),3);
    
    % Beta angle between u and w
    Segment(i).beta = mean(acosd(dot(Segment(i).Q(10:12,1,:), ...
        Segment(i).Q(1:3,1,:))),3);
    
    % Gamma angle between u and (rP - rD)
    Segment(i).gamma = mean(acosd(dot(Segment(i).Q(1:3,1,:), ...
        Segment(i).Q(4:6,1,:) - Segment(i).Q(7:9,1,:))./...
        sqrt(sum((Segment(i).Q(4:6,1,:) - ...
        Segment(i).Q(7:9,1,:)).^2))),3);
    
    % Matrix B from SCS to NSCS
    Segment(i).B = [1, ...
        Segment(i).L*cosd(Segment(i).gamma), ...
        cosd(Segment(i).beta); ...
        0, ...
        Segment(i).L*sind(Segment(i).gamma), ...
        (cosd(Segment(i).alpha) - cosd(Segment(i).beta)*cosd(Segment(i).gamma))/sind(Segment(i).gamma); ...
        0, ...
        0, ...
        sqrt(1 - cosd(Segment(i).beta)^2 - ((cosd(Segment(i).alpha) - cosd(Segment(i).beta)*cosd(Segment(i).gamma))/sind(Segment(i).gamma))^2)];
    
    % Mean coordinates of markers in (u, rP-rD, w)
    for j = 1:size(Segment(i).rM,2)
        % Projection in a non orthonormal coordinate system
        Segment(i).nM(:,j) = mean(Vnop_array3(...
            Segment(i).rM(:,j,:) - Segment(i).Q(4:6,1,:),...
            Segment(i).Q(1:3,1,:),...
            Segment(i).Q(4:6,1,:) - Segment(i).Q(7:9,1,:),...
            Segment(i).Q(10:12,1,:)),3);
    end
    
end

% -------------------------------------------------------------------------
% Universal ankle joint
% Delp (1990)
% -------------------------------------------------------------------------

% Virtual markers
% Subtalar axis direction in foot segment
ns12 = [0.781;0.600;-0.120]/norm([0.781;0.600;-0.120]); % In Delp's thesis
% Subtalar axis position in foot segment
rVs12 = [0; 0.3907 - 0.4260; 0]*(Segment(3).L/0.3907); % Scaled by shank length
% 0.3907 for shank length (gait_scale R_SHANK in SIMM)
% -0.426000 from knee joint centre to talus centre (tibia_r talus_r in SIMM)
% No rVs22
% No rVs32

% Ankle axis direction in shank segment
ns13 = [-0.105; -0.174; 0.979]/norm([-0.105; -0.174; 0.979]); % In Delp's thesis
% Ankle axis position in shank segment
% Witt regards to distal endpoints
% Ankle and subtalar axes superimposed (as well as Y-axes of foot and shank in neutral position)
rVs13 = rVs12;
% No rVs23
% No rVs33

% Coordinates of virtual markers in segment 2
% Expressed in (u2, rP2-rD2, w2)
Segment(2).nV(:,1) = inv(Segment(2).B)*rVs12; % V12: virtual marker 1 of segment 2
Segment(2).nn(:,1) = inv(Segment(2).B)*ns12;

% Coordinates of virtual markers in segment 3
% Expressed in (u3, rP3-rD3, w3)
Segment(3).nV(:,1) = [0;-1;0] + inv(Segment(3).B)*rVs13; % V13: virtual marker 1 of segment 3
Segment(3).nn(:,1) = inv(Segment(3).B)*ns13;

% Ankle complex mean angle
Nn12 = [Segment(2).nn(1,1)*eye(3),...
    (Segment(2).nn(2,1))*eye(3), ...
    - Segment(2).nn(2,1)*eye(3), ...
    Segment(2).nn(3,1)*eye(3)];
Nn13 = [Segment(3).nn(1,1)*eye(3),...
    (Segment(3).nn(2,1))*eye(3), ...
    - Segment(3).nn(2,1)*eye(3), ...
    Segment(3).nn(3,1)*eye(3)];
thetaA = mean(acos(dot(Mprod_array3(repmat(Nn13,[1,1,n]),Segment(3).Q), ...
    Mprod_array3(repmat(Nn12,[1,1,n]),Segment(2).Q))),3)*180/pi;

% -------------------------------------------------------------------------
% Parallel tibio-femoral mechanism: joint 3
% -------------------------------------------------------------------------
% Generic model of the right knee
% J.D. Feikes et al. J Biomech 36 (2003) 125-129
% In tibia SCS = femur SCS at full extension
% Origin at midpoint between medial and lateral centers
% X = -X ISB, Y = Z ISB, Z = Y ISB
% SCS adjustement: frontal plane including the 2 condyle centers
% -------------------------------------------------------------------------

% Virtual markers
rVs43 = [-5.0; -20.2; -22.3]/1000; % Medial_contact_tibia (in m)
rVs53 = [5.0; -20.2; 22.3]/1000; % Lateral_contact_tibia (in m)
rVs63 = [-0.0; -15.2; 0.0]/1000; % ACL_tibia (in m)
rVs73 = [-30.3; -23.7; -2.4]/1000; % PCL_tibia (in m)
rVs83 = [1.2; -67.2; -3.6]/1000; % MCL_tibia (in m)
rVs15 = [-5.0; 0.6; -22.3]/1000; % Medial_centre_femur (in m)
rVs25 = [5.0; -0.6; 22.3]/1000; % Lateral_centre_femur (in m)
rVs35 = [-27.8; 12.7; 5.0]/1000; % ACL_femur (in m)
rVs45 = [-20.6; -4.3; -15.7]/1000; % PCL_femur (in m)
rVs55 = [-9.7; 10.2; -42.3]/1000; % MCL_femur (in m)

% Normals (associated to segment 3)
ns23 = [0;1;0]; % Medial tibial plateau orthogonal to Y axis of femur SCS
ns33 = [0;1;0]; % Lateral tibial plateau orthogonal to Y axis of femur SCS
% Indexes incremented (ankle axis is #1)

% Spheres radii and ligament lengths: distances (associated with joint 3)
Joint(3).d(1,1) = norm(rVs15 - rVs43); % Medial sphere radius
Joint(3).d(1,2) = norm(rVs25 - rVs53); % Lateral sphere radius
Joint(3).d(1,3) = norm(rVs35 - rVs63); % ACL
Joint(3).d(1,4) = norm(rVs45 - rVs73); % PCL
Joint(3).d(1,5) = norm(rVs55 - rVs83); % MCL

% Coordinates of virtual markers in segment 3
% Expressed in (u3, rP3-rD3, w3)
Segment(3).nV(:,4) = inv(Segment(3).B)*rVs43; % V43: virtual marker 4 of segment 3
Segment(3).nV(:,5) = inv(Segment(3).B)*rVs53; % V53: virtual marker 5 of segment 3
Segment(3).nV(:,6) = inv(Segment(3).B)*rVs63; % V63 virtual marker 6 of segment 3
Segment(3).nV(:,7) = inv(Segment(3).B)*rVs73; % V73 virtual marker 7 of segment 3
Segment(3).nV(:,8) = inv(Segment(3).B)*rVs83; % V83 virtual marker 8 of segment 3

% Tibial plateau normals (associated with segment 3)
% Expressed in (u3, rP3-rD3, w3)
Segment(3).nn(:,2) = inv(Segment(3).B)*ns23;
Segment(3).nn(:,3) = inv(Segment(3).B)*ns33;
% Indexes incremented (ankle axis is #1)

% Coordinates of virtual markers in segment 5
% Expressed in (u5, rP5-rD5, w5)
Segment(5).nV(:,1) = [0;-1;0] + inv(Segment(5).B)*rVs15; % V15: virtual marker 1 of segment 5
Segment(5).nV(:,2) = [0;-1;0] + inv(Segment(5).B)*rVs25; % V25: virtual marker 2 of segment 5
Segment(5).nV(:,3) = [0;-1;0] + inv(Segment(5).B)*rVs35; % V35: virtual marker 3 of segment 5
Segment(5).nV(:,4) = [0;-1;0] + inv(Segment(5).B)*rVs45; % V45: virtual marker 4 of segment 5
Segment(5).nV(:,5) = [0;-1;0] + inv(Segment(5).B)*rVs55; % V55: virtual marker 5 of segment 5

% -------------------------------------------------------------------------
% Paralell patello-femoral mechanism: joint 4
% -------------------------------------------------------------------------
% Model of the patella hinge
% Sancisi and Parenti-Castelli
% Journal of mechanisms and robotics 3 (2011)
% In patella SCS and tibia/femur SCS at full extension
% (Model parameters before inverse fitting on cadaver kinematics)
% -------------------------------------------------------------------------

% Virtual markers
rVs93 = [20.0000; -30.0000; 5.0000]./1000 + ... % PT_tibia
    [-0.7385; -19.5520; 3.2718]./1000; % Translation from femur SCS to tibia SCS (values provided by Sancisi)
rVs14 = [-42.48; 6.18; 0]./1000; % P_axis_patella
% rVs24 = rD4:   PT_patella
rVs65 = [4.77; 10.47; -2.68]./1000;  % P_axis_femur

% Axis orientation
delta_patella = -0.13; % Azimuth of the hinge axis, expressed in the patella SCS
nu_patella = 0.19; % Altitude of the hinge axis, expressed in the patella SCS
ns14 = [sin(nu_patella); ... % n_axis_patella : orientation of the hinge axis, expressed in the patella SCS
    cos(nu_patella)*sin(delta_patella); ...
    cos(nu_patella)*cos(delta_patella)];
ns14 = ns14/norm(ns14); % Unitary vector
delta_femur = 0.09; % Azimuth of the hinge axis, expressed in the femur SCS
nu_femur = 0.09; % Altitude of the hinge axis, expressed in the femur SCS
ns15 = [sin(nu_femur); ...  % n_axis_femur : orientation of the hinge axis, expressed in the femur SCS
    cos(nu_femur)*sin(delta_femur); ...
    cos(nu_femur)*cos(delta_femur)];
ns15 = ns15/norm(ns15); % Unitary vector

% Orientation of hinge axis
thetaP1 = acosd(dot(ns15,[1.0, 0.0, 0.0])); % Angle between n_axis_femur and X axis of femur SCS
thetaP2 = acosd(dot(ns15,[0.0, 1.0, 0.0])); % Angle between n_axis_femur and Y axis of femur SCS

% Patellar tendon length: distances associated with joint 4
Joint(4).d(1,1) = 48.67/1000;

% Coordinates of virtual markers in segment 3
% Expressed in (u3, rP3-rD3, w3)
Segment(3).nV(:,9) = inv(Segment(3).B)*rVs93; % V93: virtual marker 9 of segment 3

% Coordinates of virtual markers in segment 4
% Expressed in (u4, rP4-rD4, w4)
Segment(4).nV(:,1) = inv(Segment(4).B)*rVs14; % V14: virtual marker 1 of segment 4
% Orientation (associated with segment 4)
Segment(4).nn(:,1) = inv(Segment(4).B)*ns14; % n14: virtual normal 1 of segment 4

% Coordinates of virtual markers in segment 5
% Expressed in (u5, rP5-rD5, w5)
Segment(5).nV(:,6) = [0;-1;0] + inv(Segment(5).B)*rVs65; % V65: virtual marker 6 of segment 5
% Orientation (associated with segment 5)
Segment(5).nn(:,1) = inv(Segment(5).B)*ns15; % n15: virtual normal 1 of segment 5

% -------------------------------------------------------------------------
% Spherical hip joint: joint 5
% -------------------------------------------------------------------------
% Hip virtual marker mean coordinates (rV1 = rP6)
% Expressed in  in (u6, rP6-rD6, w6)
Segment(6).nV(:,1) = mean(Vnop_array3(...
    Segment(5).Q(4:6,1,:) - Segment(6).Q(4:6,1,:),...
    Segment(6).Q(1:3,1,:),...
    Segment(6).Q(4:6,1,:) - Segment(6).Q(7:9,1,:),...
    Segment(6).Q(10:12,1,:)),3);

% -------------------------------------------------------------------------
% Interpolation matrices
% -------------------------------------------------------------------------
% Segment 2
NV12 = [Segment(2).nV(1,1)*eye(3),...
    (1 + Segment(2).nV(2,1))*eye(3), ...
    - Segment(2).nV(2,1)*eye(3), ...
    Segment(2).nV(3,1)*eye(3)];
Nn12 = [Segment(2).nn(1,1)*eye(3),...
    (Segment(2).nn(2,1))*eye(3), ...
    - Segment(2).nn(2,1)*eye(3), ...
    Segment(2).nn(3,1)*eye(3)];

% Segment 3
NV13 = [Segment(3).nV(1,1)*eye(3),...
    (1 + Segment(3).nV(2,1))*eye(3), ...
    - Segment(3).nV(2,1)*eye(3), ...
    Segment(3).nV(3,1)*eye(3)];
Nn13 = [Segment(3).nn(1,1)*eye(3),...
    (Segment(3).nn(2,1))*eye(3), ...
    - Segment(3).nn(2,1)*eye(3), ...
    Segment(3).nn(3,1)*eye(3)];

% NV23 = [Segment(3).nV(1,2)*eye(3),...
%     (1 + Segment(3).nV(2,2))*eye(3), ...
%     - Segment(3).nV(2,2)*eye(3), ...
%     Segment(3).nV(3,2)*eye(3)];
% NV33 = [Segment(3).nV(1,3)*eye(3),...
%     (1 + Segment(3).nV(2,3))*eye(3), ...
%     - Segment(3).nV(2,3)*eye(3), ...
%     Segment(3).nV(3,3)*eye(3)];
NV43 = [Segment(3).nV(1,4)*eye(3),...
    (1 + Segment(3).nV(2,4))*eye(3), ...
    - Segment(3).nV(2,4)*eye(3), ...
    Segment(3).nV(3,4)*eye(3)];
NV53 = [Segment(3).nV(1,5)*eye(3),...
    (1 + Segment(3).nV(2,5))*eye(3), ...
    - Segment(3).nV(2,5)*eye(3), ...
    Segment(3).nV(3,5)*eye(3)];
NV63 = [Segment(3).nV(1,6)*eye(3),...
    (1 + Segment(3).nV(2,6))*eye(3), ...
    - Segment(3).nV(2,6)*eye(3), ...
    Segment(3).nV(3,6)*eye(3)];
NV73 = [Segment(3).nV(1,7)*eye(3),...
    (1 + Segment(3).nV(2,7))*eye(3), ...
    - Segment(3).nV(2,7)*eye(3), ...
    Segment(3).nV(3,7)*eye(3)];
NV83 = [Segment(3).nV(1,8)*eye(3),...
    (1 + Segment(3).nV(2,8))*eye(3), ...
    - Segment(3).nV(2,8)*eye(3), ...
    Segment(3).nV(3,8)*eye(3)];
NV93 = [Segment(3).nV(1,9)*eye(3),...
    (1 + Segment(3).nV(2,9))*eye(3), ...
    - Segment(3).nV(2,9)*eye(3), ...
    Segment(3).nV(3,9)*eye(3)];
Nn23 = [Segment(3).nn(1,2)*eye(3),...
    (Segment(3).nn(2,2))*eye(3), ...
    - Segment(3).nn(2,2)*eye(3), ...
    Segment(3).nn(3,2)*eye(3)];
Nn33 = [Segment(3).nn(1,3)*eye(3),...
    (Segment(3).nn(2,3))*eye(3), ...
    - Segment(3).nn(2,3)*eye(3), ...
    Segment(3).nn(3,3)*eye(3)];

% Segment 4
NV14 = [Segment(4).nV(1,1)*eye(3),...
    (1 + Segment(4).nV(2,1))*eye(3), ...
    - Segment(4).nV(2,1)*eye(3), ...
    Segment(4).nV(3,1)*eye(3)];
Nn14 = [Segment(4).nn(1,1)*eye(3),...
    (Segment(4).nn(2,1))*eye(3), ...
    - Segment(4).nn(2,1)*eye(3), ...
    Segment(4).nn(3,1)*eye(3)];

% Segment 5
NV15 = [Segment(5).nV(1,1)*eye(3),...
    (1 + Segment(5).nV(2,1))*eye(3), ...
    - Segment(5).nV(2,1)*eye(3), ...
    Segment(5).nV(3,1)*eye(3)];
NV25 = [Segment(5).nV(1,2)*eye(3),...
    (1 + Segment(5).nV(2,2))*eye(3), ...
    - Segment(5).nV(2,2)*eye(3), ...
    Segment(5).nV(3,2)*eye(3)];
NV35 = [Segment(5).nV(1,3)*eye(3),...
    (1 + Segment(5).nV(2,3))*eye(3), ...
    - Segment(5).nV(2,3)*eye(3), ...
    Segment(5).nV(3,3)*eye(3)];
NV45 = [Segment(5).nV(1,4)*eye(3),...
    (1 + Segment(5).nV(2,4))*eye(3), ...
    - Segment(5).nV(2,4)*eye(3), ...
    Segment(5).nV(3,4)*eye(3)];
NV55 = [Segment(5).nV(1,5)*eye(3),...
    (1 + Segment(5).nV(2,5))*eye(3), ...
    - Segment(5).nV(2,5)*eye(3), ...
    Segment(5).nV(3,5)*eye(3)];
NV65 = [Segment(5).nV(1,6)*eye(3),...
    (1 + Segment(5).nV(2,6))*eye(3), ...
    - Segment(5).nV(2,6)*eye(3), ...
    Segment(5).nV(3,6)*eye(3)];

% Segment 6
NV16 = [Segment(6).nV(1,1)*eye(3),...
    (1 + Segment(6).nV(2,1))*eye(3), ...
    - Segment(6).nV(2,1)*eye(3), ...
    Segment(6).nV(3,1)*eye(3)];


%% -------------------------------------------------------------------------
% Run optimisation
% -------------------------------------------------------------------------
% Initial guess for Lagrange multipliers
lambdar = zeros(30,1,n); % 5 segments x 6 constraints per segment
lambdakA = zeros(4,1,n);
lambdakK = zeros(5,1,n);
lambdakP = zeros(6,1,n);
lambdakH = zeros(3,1,n);

% Initial value of the objective function
F = 1;
% Iteration number
step = 0;

% Newton-Raphson
while max(permute(sqrt(sum(F.^2)),[3,2,1])) > 10e-12 && step < 20
    
    % Iteration number
    step = step + 1 % Display
    
    % Initialisation
    phik = []; % Vector of kinematic constraints
    Kk = [];  % Jacobian of kinematic constraints
    phir = []; % Vector of rigid body constraints
    Kr = []; % Jacobian of rigid body constraints
    dKlambdardQ = []; % Partial derivative of Jacobian * Lagrange multipliers
    phim = []; % Vector of driving constraints
    Km = []; % Jacobian of driving constraints
    
    % Universal ankle joint
    % ---------------------------------------------------------------------
    
    % Position of virtual markers
    rV12 = Mprod_array3(repmat(NV12,[1,1,n]),Segment(2).Q);
    n12 = Mprod_array3(repmat(Nn12,[1,1,n]),Segment(2).Q);  % Subtalar axis
    rV13 = Mprod_array3(repmat(NV13,[1,1,n]),Segment(3).Q);
    n13 = Mprod_array3(repmat(Nn13,[1,1,n]),Segment(3).Q);  % Ankle axis
    
    % Vector of kinematic constraints
    % rV13 - rV12 = 0 (Lateral_centre_calcaneus = Lateral_contact_tibia)
    % n13.n12 - cos(thetaA) = 0
    phikA = [rV13 - rV12;...
        dot(n13,n12) - repmat(cosd(thetaA),[1,1,n])];
    
    % Jacobian of kinematic constraints
    KkA = zeros(4,5*12,n); % Initialisation
    KkA(1:3,1:12,:) = - repmat(NV12,[1,1,n]);
    KkA(1:3,13:24,:) = repmat(NV13,[1,1,n]);
    KkA(4,1:12,:) = Mprod_array3(permute(n13,[2,1,3]),repmat(Nn12,[1,1,n]));
    KkA(4,13:24,:) = Mprod_array3(permute(n12,[2,1,3]),repmat(Nn13,[1,1,n]));
    % with transpose = permute( ,[2,1,3])
    
    % Joint structure
    Joint(2).Kk = KkA;
    
    % Partial derivative of Jacobian * Lagrange multipliers
    dKlambdakAdQ = zeros(5*12,5*12,n); % Initialisation
    dKlambdakAdQ(1:12,13:24,:) = Mprod_array3(lambdakA(4,1,:),repmat(Nn12'*Nn13,[1,1,n]));
    dKlambdakAdQ(13:24,1:12,:) = Mprod_array3(lambdakA(4,1,:),repmat(Nn13'*Nn12,[1,1,n]));


    % Parallel tibio-femoral mechanism
    % ---------------------------------------------------------------------
    % Position of virtual markers
    rV43 = Mprod_array3(repmat(NV43,[1,1,n]),Segment(3).Q); % Medial_contact_tibia
    rV53 = Mprod_array3(repmat(NV53,[1,1,n]),Segment(3).Q); % Lateral_contact_tibia
    rV63 = Mprod_array3(repmat(NV63,[1,1,n]),Segment(3).Q); % ACL_tibia
    rV73 = Mprod_array3(repmat(NV73,[1,1,n]),Segment(3).Q); % PCL_tibia
    rV83 = Mprod_array3(repmat(NV83,[1,1,n]),Segment(3).Q); % MCL_tibia
    rV15 = Mprod_array3(repmat(NV15,[1,1,n]),Segment(5).Q); % Medial_centre_femur
    rV25 = Mprod_array3(repmat(NV25,[1,1,n]),Segment(5).Q); % Lateral_centre_femur
    rV35 = Mprod_array3(repmat(NV35,[1,1,n]),Segment(5).Q); % ACL_femur
    rV45 = Mprod_array3(repmat(NV45,[1,1,n]),Segment(5).Q); % PCL_femur
    rV55 = Mprod_array3(repmat(NV55,[1,1,n]),Segment(5).Q); % MCL_femur
    
    % Direction of normals
    n23 = Mprod_array3(repmat(Nn23,[1,1,n]),Segment(3).Q); % Medial_tibia
    n33 = Mprod_array3(repmat(Nn33,[1,1,n]),Segment(3).Q); % Lateral_tibia
    % Indexes incremented (ankle axis is #1)
    
    % Vector of kinematic constraints
    phikK = [dot((rV15 - rV43),n23) - repmat(Joint(3).d(1,1),[1,1,n]);...
        dot((rV25 - rV53),n33) - repmat(Joint(3).d(1,2),[1,1,n]);...
        dot((rV35 - rV63),(rV35 - rV63)) - repmat((Joint(3).d(1,3))^2,[1,1,n]);...
        dot((rV45 - rV73),(rV45 - rV73)) - repmat((Joint(3).d(1,4))^2,[1,1,n]);...
        dot((rV55 - rV83),(rV55 - rV83)) - repmat((Joint(3).d(1,5))^2,[1,1,n])];
    
    % Jacobian of kinematic constraints
    KkK = zeros(5,5*12,n);
    KkK(1,13:24,:) = - Mprod_array3(permute(n23,[2,1,3]),repmat(NV43,[1,1,n])) ...
        + Mprod_array3(permute(rV15 - rV43,[2,1,3]),repmat(Nn23,[1,1,n]));
    KkK(1,37:48,:) = Mprod_array3(permute(n23,[2,1,3]),repmat(NV15,[1,1,n]));
    KkK(2,13:24,:) = - Mprod_array3(permute(n33,[2,1,3]),repmat(NV53,[1,1,n])) ...
        + Mprod_array3(permute(rV25 - rV53,[2,1,3]),repmat(Nn33,[1,1,n]));
    KkK(2,37:48,:) = Mprod_array3(permute(n33,[2,1,3]),repmat(NV25,[1,1,n]));
    KkK(3,13:24,:) = - 2*Mprod_array3(permute(rV35 - rV63,[2,1,3]),repmat(NV63,[1,1,n]));
    KkK(3,37:48,:) = 2*Mprod_array3(permute(rV35 - rV63,[2,1,3]),repmat(NV35,[1,1,n]));
    KkK(4,13:24,:) = - 2*Mprod_array3(permute(rV45 - rV73,[2,1,3]),repmat(NV73,[1,1,n]));
    KkK(4,37:48,:) = 2*Mprod_array3(permute(rV45 - rV73,[2,1,3]),repmat(NV45,[1,1,n]));
    KkK(5,13:24,:) = - 2*Mprod_array3(permute(rV55 - rV83,[2,1,3]),repmat(NV83,[1,1,n]));
    KkK(5,37:48,:) = 2*Mprod_array3(permute(rV55 - rV83,[2,1,3]),repmat(NV55,[1,1,n]));
    % with transpose = permute( ,[2,1,3])
    % Joint structure
    Joint(3).Kk = KkK;
    
    % Partial derivative of Jacobian * Lagrange multipliers
    dKlambdakKdQ = zeros(5*12,5*12,n); % Initialisation
    dKlambdakKdQ(13:24,13:24,:) = - Mprod_array3(lambdakK(1,1,:),repmat(NV43'*Nn23 + Nn23'*NV43,[1,1,n])) ...
        - Mprod_array3(lambdakK(2,1,:),repmat(NV53'*Nn33 + Nn33'*NV53,[1,1,n])) ...
        + Mprod_array3(lambdakK(3,1,:),repmat(2*NV63'*NV63,[1,1,n])) ...
        + Mprod_array3(lambdakK(4,1,:),repmat(2*NV73'*NV73,[1,1,n])) ...
        + Mprod_array3(lambdakK(5,1,:),repmat(2*NV83'*NV83,[1,1,n]));
    dKlambdakKdQ(13:24,37:48,:) = Mprod_array3(lambdakK(1,1,:),repmat(Nn23'*NV15,[1,1,n])) ...
        + Mprod_array3(lambdakK(2,1,:),repmat(Nn33'*NV25,[1,1,n])) ...
        - Mprod_array3(lambdakK(3,1,:),repmat(2*NV63'*NV35,[1,1,n])) ...
        - Mprod_array3(lambdakK(4,1,:),repmat(2*NV73'*NV45,[1,1,n])) ...
        - Mprod_array3(lambdakK(5,1,:),repmat(2*NV83'*NV55,[1,1,n]));
    dKlambdakKdQ(37:48,13:24,:) = permute(dKlambdakKdQ(13:24,37:48,:),[2,1,3]); % Symetrical
    % with transpose = permute( ,[2,1,3])
    dKlambdakKdQ(37:48,37:48,:) = Mprod_array3(lambdakK(3,1,:),repmat(2*NV35'*NV35,[1,1,n])) ...
        + Mprod_array3(lambdakK(4,1,:),repmat(2*NV45'*NV45,[1,1,n])) ...
        + Mprod_array3(lambdakK(5,1,:),repmat(2*NV55'*NV55,[1,1,n]));
    
    % Parallel patello-femoral mechanism
    % ---------------------------------------------------------------------
    % Position of virtual markers
    rV93 = Mprod_array3(repmat(NV93,[1,1,n]),Segment(3).Q); % PT_tibia
    rV14 = Mprod_array3(repmat(NV14,[1,1,n]),Segment(4).Q); % P_axis_patella
    rV65 = Mprod_array3(repmat(NV65,[1,1,n]),Segment(5).Q); % P_axis_femur
    
    % Direction of normals
    n14 = Mprod_array3(repmat(Nn14,[1,1,n]),Segment(4).Q); % n_axis_patella
    %     n15 = Mprod_array3(repmat(Nn15,[1,1,n]),Segment(5).Q);
    
    % Vector of kinematic constraints
    phikP = [rV65 - rV14; ...
        dot(Segment(5).Q(1:3,1,:),n14)- cosd(thetaP1); ...
        dot((Segment(5).Q(4:6,1,:)- Segment(5).Q(7:9,1,:)), n14) - Segment(5).L*cosd(thetaP2); ...
        dot((Segment(4).Q(7:9,1,:) - rV93),(Segment(4).Q(7:9,1,:) - rV93)) - (repmat((Joint(4).d(1,1))^2,[1,1,n]))];
    
    % Jacobian of kinematic constraints
    KkP = zeros(6,5*12,n); % Initialisation
    KkP(1:3,25:36,:) = repmat(-NV14,[1,1,n]);
    KkP(1:3,37:48,:) = repmat(NV65,[1,1,n]);
    KkP(4,25:36,:) = Mprod_array3(permute(Segment(5).Q(1:3,1,:),[2,1,3]),repmat(Nn14,[1,1,n]));
    KkP(4,37:48,:) = [permute(n14,[2,1,3]),zeros(1,3,n),zeros(1,3,n),zeros(1,3,n)];
    KkP(5,25:36,:) = Mprod_array3(permute(Segment(5).Q(4:6,1,:) - Segment(5).Q(7:9,1,:),[2,1,3]),repmat(Nn14,[1,1,n]));
    KkP(5,37:48,:) = [zeros(1,3,n),permute(n14,[2,1,3]),permute(-n14,[2,1,3]),zeros(1,3,n)];
    KkP(6,13:24,:) = - Mprod_array3(permute(Segment(4).Q(7:9,1,:) - rV93,[2,1,3]),repmat(2*NV93,[1,1,n]));
    KkP(6,31:33,:) = Mprod_array3(permute(Segment(4).Q(7:9,1,:) - rV93,[2,1,3]),repmat(2*eye(3),[1,1,n]));
    % with transpose = permute( ,[2,1,3])
    % Joint structure
    Joint(4).Kk = KkP;
    
    % Partial derivative of Jacobian * Lagrange multipliers
    dKlambdakPdQ = zeros(5*12,5*12,n); % Initialisation
    dKlambdakPdQ(13:24,13:24,:) = Mprod_array3(lambdakP(6,1,:),repmat(2*NV93'*NV93,[1,1,n]));
    dKlambdakPdQ(19:21,25:36,:) = - Mprod_array3(lambdakP(6,1,:),repmat(2*eye(3)*NV93,[1,1,n]));
    dKlambdakPdQ(31:33,13:24,:) = - Mprod_array3(lambdakP(6,1,:),repmat(2*eye(3)*NV93,[1,1,n]));
    dKlambdakPdQ(31:33,31:33,:) = Mprod_array3(lambdakP(6,1,:),repmat(-2*eye(3),[1,1,n]));
    dKlambdakPdQ(25:36,37:48,:) = [Mprod_array3(lambdakP(4,1,:),repmat(Nn14,[1,1,n]));...
        Mprod_array3(lambdakP(5,1,:),repmat(Nn14,[1,1,n]));...
        Mprod_array3(lambdakP(5,1,:),repmat(-Nn14,[1,1,n]));...
        repmat(zeros(3,12),[1,1,n])];
    dKlambdakPdQ(37:48,25:36,:) = [Mprod_array3(lambdakP(4,1,:),repmat(Nn14,[1,1,n]));...
        Mprod_array3(lambdakP(5,1,:),repmat(Nn14,[1,1,n]));...
        Mprod_array3(lambdakP(5,1,:),repmat(-Nn14,[1,1,n]));...
        repmat(zeros(3,12),[1,1,n])];
    
    % Spherical hip joint
    % ---------------------------------------------------------------------
    
    % Position of virtual markers
    rV16 = Mprod_array3(repmat(NV16,[1,1,n]),Segment(6).Q); % HJC_pelvis
    
    % Vector of kinematic constraints
    phikH = rV16 - Segment(5).Q(4:6,1,:);
    
    % Jacobian of kinematic constraints
    KkH = zeros(3,5*12,n); % Initialisation
    KkH(1:3,40:42,:) = repmat(-eye(3),[1,1,n]);
    KkH(1:3,49:60,:) = repmat(NV16,[1,1,n]);
    % Joint structure
    Joint(5).Kk = KkH;
    
    % Partial derivative of Jacobian * Lagrange multipliers
    dKlambdakHdQ = zeros(5*12,5*12,n); % Initialisation
    
    % Assembly
    % ---------------------------------------------------------------------
    phik = [phikA;phikK;phikP;phikH];
    Kk = [KkA;KkK;KkP;KkH];
    
    % Rigid body constraints
    % ---------------------------------------------------------------------
    for i = 2:6 % From i = 2 (Foot) to i = 6 (Pelvis)
        
        % Vector of rigid body constraints
        ui = Segment(i).Q(1:3,1,:);
        vi = Segment(i).Q(4:6,1,:) - Segment(i).Q(7:9,1,:);
        wi = Segment(i).Q(10:12,1,:);
        phiri = [dot(ui,ui) - ones(1,1,n);...
            dot(ui,vi) - repmat(Segment(i).L*cosd(Segment(i).gamma),[1,1,n]); ...
            dot(ui,wi) - repmat(cosd(Segment(i).beta),[1,1,n]); ...
            dot(vi,vi) - repmat(Segment(i).L^2,[1,1,n]);
            dot(vi,wi) - repmat(Segment(i).L*cosd(Segment(i).alpha),[1,1,n]);
            dot(wi,wi) - ones(1,1,n)];
        
        % Jacobian of rigid body constraints
        Kri = zeros(6,5*12,n); % Initialisation
        Kri(1:6,(i-2)*12+1:(i-2)*12+12,:) = permute(...
            [    2*ui,       vi,           wi,     zeros(3,1,n),zeros(3,1,n),zeros(3,1,n); ...
            zeros(3,1,n),    ui,      zeros(3,1,n),    2*vi,         wi,     zeros(3,1,n); ...
            zeros(3,1,n),   -ui,      zeros(3,1,n),   -2*vi,        -wi,     zeros(3,1,n); ...
            zeros(3,1,n),zeros(3,1,n),     ui,     zeros(3,1,n),     vi,         2*wi],[2,1,3]);
        % with transpose = permute( ,[2,1,3])
        % Segment structure
        Segment(i).Kr = Kri;
        
        % Partial derivative of Jacobian * Lagrange multipliers
        dKlambdaridQ = zeros(12,5*12,n); % Initialisation
        lambdari = lambdar((i-2)*6+1:(i-2)*6+6,1,:); % Extraction
        dKlambdaridQ(1:12,(i-2)*12+1:(i-2)*12+12,:) = ...
            [Mprod_array3(lambdari(1,1,:),repmat(2*eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(2,1,:),repmat(eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(2,1,:),repmat(-1*eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(3,1,:),repmat(eye(3),[1,1,n])); ...
            Mprod_array3(lambdari(2,1,:),repmat(eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(4,1,:),repmat(2*eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(4,1,:),repmat(-2*eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(5,1,:),repmat(eye(3),[1,1,n])); ...
            Mprod_array3(lambdari(2,1,:),repmat(-1*eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(4,1,:),repmat(-2*eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(4,1,:),repmat(2*eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(5,1,:),repmat(-1*eye(3),[1,1,n])); ...
            Mprod_array3(lambdari(3,1,:),repmat(eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(5,1,:),repmat(eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(5,1,:),repmat(-1*eye(3),[1,1,n])), ...
            Mprod_array3(lambdari(6,1,:),repmat(2*eye(3),[1,1,n]))];
        
        % Vector and Jacobian of driving constraints
        Kmi = zeros(size(Segment(i).rM,2)*3,5*12,n); % Initialisation
        phimi = []; % Initialisation
        
        for j = 1:size(Segment(i).rM,2)
            % Interpolation matrix
            NMij = [Segment(i).nM(1,j)*eye(3),...
                (1 + Segment(i).nM(2,j))*eye(3), ...
                - Segment(i).nM(2,j)*eye(3), ...
                Segment(i).nM(3,j)*eye(3)];
            % Vector of driving constraints
            phimi((j-1)*3+1:(j-1)*3+3,1,:) = Segment(i).rM(:,j,:) ...
                - Mprod_array3(repmat(NMij,[1,1,n]),Segment(i).Q);
            % Jacobian of driving contraints
            Kmi((j-1)*3+1:(j-1)*3+3,(i-2)*12+1:(i-2)*12+12,:) = ...
                repmat(-NMij,[1,1,n]);
        end
        
        % Assembly
        phir = [phir;phiri];
        Kr = [Kr;Kri];
        dKlambdardQ = [dKlambdardQ;dKlambdaridQ];
        phim = [phim;phimi];
        Km = [Km;Kmi];
        
    end
    
    % Display errors
    % ---------------------------------------------------------------------
    Mean_phik = mean(Mprod_array3(permute(phik,[2,1,3]),phik),3)
    Mean_phir = mean(Mprod_array3(permute(phir,[2,1,3]),phir),3)
    Mean_phim = mean(Mprod_array3(permute(phim,[2,1,3]),phim),3)
    
    % Compute dX
    % ---------------------------------------------------------------------
    % dX = inv(-dF/dX)*F(X)
    % F(X) = [Km'*phim + [Kk;Kr]'*[lambdak;lambdar];[phik;phir]]
    % X = [Q;[lambdak;lambdar]]
    F = [Mprod_array3(permute(Km,[2,1,3]),phim) ...
        + Mprod_array3(permute([Kk;Kr],[2,1,3]), [lambdakA;lambdakK;lambdakP;lambdakH;lambdar]); ...
        [phik;phir]];
    % with transpose = permute( ,[2,1,3])
    dKlambdadQ = dKlambdakAdQ + dKlambdakKdQ + dKlambdakPdQ + dKlambdakHdQ + dKlambdardQ;
    dFdX = [Mprod_array3(permute(Km,[2,1,3]),Km) + dKlambdadQ, permute([Kk;Kr],[2,1,3]); ...
        [Kk;Kr],zeros(size([Kk;Kr],1),size([Kk;Kr],1),n)];
    dX = Mprod_array3(Minv_array3(-dFdX),F);
    
    % Extraction from X
    % ---------------------------------------------------------------------
    Segment(2).Q = Segment(2).Q + dX(1:12,1,:);     % Foot
    Segment(3).Q = Segment(3).Q + dX(13:24,1,:);    % Shank
    Segment(4).Q = Segment(4).Q + dX(25:36,1,:);    % Patella
    Segment(5).Q = Segment(5).Q + dX(37:48,1,:);    % Thigh
    Segment(6).Q = Segment(6).Q + dX(49:60,1,:);    % Pelvis
    lambdakA = lambdakA + dX(61:64,1,:);            % Ankle
    lambdakK = lambdakK + dX(65:69,1,:);            % Knee
    lambdakP = lambdakP + dX(70:75,1,:);            % Patella
    lambdakH = lambdakH + dX(76:78,1,:);            % Hip
    lambdar = lambdar + dX(79:end,1,:);             % Rigid constraints
    
end


%% -------------------------------------------------------------------------
% Inconsistent velocities and accelerations
% -------------------------------------------------------------------------

% Q
Q = [Segment(2).Q; ...  % Foot
    Segment(3).Q; ...   % Shank
    Segment(4).Q; ...   % Patella
    Segment(5).Q; ...   % Thigh
    Segment(6).Q];      % Pelvis

% dQdt with filtering
dQdt = Vfilt_array3(Derive_array3(Q,1/f),f,fc);
% Extraction from dQdt
Segment(2).dQdt = dQdt(1:12,1,:); % Foot
Segment(3).dQdt = dQdt(13:24,1,:); % Shank
Segment(4).dQdt = dQdt(25:36,1,:); % Patella
Segment(5).dQdt = dQdt(37:48,1,:); % Thigh
Segment(6).dQdt = dQdt(49:60,1,:); % Pelvis

% d2Qdt2 with filtering
d2Qdt2 = Vfilt_array3(Derive_array3(dQdt,1/f),f,fc);
% Extraction from d2Qdt2
Segment(2).d2Qdt2 = d2Qdt2(1:12,1,:); % Foot
Segment(3).d2Qdt2 = d2Qdt2(13:24,1,:); % Shank
Segment(4).d2Qdt2 = d2Qdt2(25:36,1,:); % Patella
Segment(5).d2Qdt2 = d2Qdt2(37:48,1,:); % Thigh
Segment(6).d2Qdt2 = d2Qdt2(49:60,1,:); % Pelvis


%% -------------------------------------------------------------------------
% Modified Jacobian of the kinematic constraint for spherical joints
% Constraints projected on the axes of the distal segment
% -------------------------------------------------------------------------

% Ankle
% X2
NX2 = [eye(3),zeros(3,9)]; % X2 = u2
X2 = Segment(2).Q(1:3,1,:); % u2 = X2 of SCS
% Y2
inB2 = inv(Segment(2).B);
NY2 = [inB2(1,2)*eye(3),inB2(2,2)*eye(3),-inB2(2,2)*eye(3),inB2(3,2)*eye(3)];
Y2 = Mprod_array3(repmat(NY2,[1,1,n]),Segment(2).Q); % Y2 of SCS
% Z2
NZ2 = [inB2(1,3)*eye(3),inB2(2,3)*eye(3),-inB2(2,3)*eye(3),inB2(3,3)*eye(3)];
Z2 = Mprod_array3(repmat(NZ2,[1,1,n]),Segment(2).Q); % Z2 of SCS
%
Joint(2).Kk(1,1:12,:) = -Mprod_array3(permute(X2,[2,1,3]),repmat(NV12,[1,1,n])) + ...
    Mprod_array3(permute(rV13 - rV12,[2,1,3]),repmat(NX2,[1,1,n]));
Joint(2).Kk(1,13:24,:) = Mprod_array3(permute(X2,[2,1,3]),repmat(NV13,[1,1,n]));
Joint(2).Kk(2,1:12,:) = -Mprod_array3(permute(Y2,[2,1,3]),repmat(NV12,[1,1,n])) + ...
    Mprod_array3(permute(rV13 - rV12,[2,1,3]),repmat(NY2,[1,1,n]));
Joint(2).Kk(2,13:24,:) = Mprod_array3(permute(Y2,[2,1,3]),repmat(NV13,[1,1,n]));
Joint(2).Kk(3,1:12,:) = -Mprod_array3(permute(Z2,[2,1,3]),repmat(NV12,[1,1,n])) + ...
    Mprod_array3(permute(rV13 - rV12,[2,1,3]),repmat(NZ2,[1,1,n]));
Joint(2).Kk(3,13:24,:) = Mprod_array3(permute(Z2,[2,1,3]),repmat(NV13,[1,1,n]));

% Patello-femoral
% X4
NX4 = [eye(3),zeros(3,9)]; % X4 = u4
X4 = Segment(4).Q(1:3,1,:); % u4 = X4 of SCS
% Y4
inB4 = inv(Segment(4).B);
NY4 = [inB4(1,2)*eye(3),inB4(2,2)*eye(3),-inB4(2,2)*eye(3),inB4(3,2)*eye(3)];
Y4 = Mprod_array3(repmat(NY4,[1,1,n]),Segment(4).Q); % Y4 of SCS
% Z4
NZ4 = [inB4(1,3)*eye(3),inB4(2,3)*eye(3),-inB4(2,3)*eye(3),inB4(3,3)*eye(3)];
Z4 = Mprod_array3(repmat(NZ4,[1,1,n]),Segment(4).Q); % Z4 of SCS
%
Joint(4).Kk(1,25:36,:) = -Mprod_array3(permute(X4,[2,1,3]),repmat(NV14,[1,1,n])) + ...
    Mprod_array3(permute(rV65 - rV14,[2,1,3]),repmat(NX4,[1,1,n]));
Joint(4).Kk(1,37:48,:) = Mprod_array3(permute(X4,[2,1,3]),repmat(NV65,[1,1,n]));
Joint(4).Kk(2,25:36,:) = -Mprod_array3(permute(Y4,[2,1,3]),repmat(NV14,[1,1,n])) + ...
    Mprod_array3(permute(rV65 - rV14,[2,1,3]),repmat(NY4,[1,1,n]));
Joint(4).Kk(2,37:48,:) = Mprod_array3(permute(Y4,[2,1,3]),repmat(NV65,[1,1,n]));
Joint(4).Kk(3,25:36,:) = -Mprod_array3(permute(Z4,[2,1,3]),repmat(NV14,[1,1,n])) + ...
    Mprod_array3(permute(rV65 - rV14,[2,1,3]),repmat(NZ4,[1,1,n]));
Joint(4).Kk(3,37:48,:) = Mprod_array3(permute(Z4,[2,1,3]),repmat(NV65,[1,1,n]));

% Hip
% X5
NX5 = [eye(3),zeros(3,9)]; % X5 = u5
X5 = Segment(4).Q(1:3,1,:); % u5 = X5 of SCS
% Y5
inB5 = inv(Segment(5).B);
NY5 = [inB5(1,2)*eye(3),inB5(2,2)*eye(3),-inB5(2,2)*eye(3),inB5(3,2)*eye(3)];
Y5 = Mprod_array3(repmat(NY5,[1,1,n]),Segment(5).Q); % Y5 of SCS
% Z5
NZ5 = [inB5(1,3)*eye(3),inB5(2,3)*eye(3),-inB5(2,3)*eye(3),inB5(3,3)*eye(3)];
Z5 = Mprod_array3(repmat(NZ5,[1,1,n]),Segment(5).Q); % Z5 of SCS
%
Joint(5).Kk(1,37:48,:) = [zeros(1,3,n),permute(-X5,[2,1,3]),zeros(1,6,n)] + ...
    Mprod_array3(permute(rV16 - Segment(5).Q(4:6,1,:),[2,1,3]),repmat(NX5,[1,1,n]));
Joint(5).Kk(2,37:48,:) = [zeros(1,3,n),permute(-Y5,[2,1,3]),zeros(1,6,n)] + ...
    Mprod_array3(permute(rV16 - Segment(5).Q(4:6,1,:),[2,1,3]),repmat(NY5,[1,1,n]));
Joint(5).Kk(3,37:48,:) = [zeros(1,3,n),permute(-Z5,[2,1,3]),zeros(1,6,n)] + ...
    Mprod_array3(permute(rV16 - Segment(5).Q(4:6,1,:),[2,1,3]),repmat(NZ5,[1,1,n]));

