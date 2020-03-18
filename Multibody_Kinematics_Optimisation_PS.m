% FUNCTION
% Multibody_Kinematics_Optimisation_PS.m
%__________________________________________________________________________
%
% PURPOSE
% Multibody optimisation of the postions of foot, shank, patella, thigh
% and computation of their consistent velocities and accelerations with
% universal ankle, tibio-femoral contact points, paralell patello-femoral 
% mechanism and spherial hip joint
%
% SYNOPSIS
% [Segment,Joint] = Multibody_Kinematics_Optimisation_GC5(Segment,Joint,f)
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
% Computation of the consistent velocities and accelerations by projection
% on the nullspace of the Jacobian and time derivative of Jacobian matrices
%
% Joint models are universal ankle, tibio-femoral contact points, paralell 
% patello-femoral mechanism, spherial hip joint
%
% REFERENCES
% F. Moissenet, L. Cheze, R. Dumas, Anatomical kinematic constraints: 
% consequences on musculo-tendon forces and joint reactions. Multibody 
% System Dynamics 2012, 28:125�141.
% S. Duprey, L. Cheze, R. Dumas. Influence of joint constraints on lower 
% limb kinematics estimation from skin markers using global optimization.
% Journal of Biomechancis 2010, 43(14):2858-62.
% A. Zeighami, R. Aissaoui, R. Dumas. Knee medial and lateral contact 
% forces in a musculoskeletal model with subject-specific contact point 
% trajectories. Journal of Biomechanics 2018, 69: 138�145.
% K.M. Varadarajan, A.L. Moynihan , D. D'Lima , C.W. Colwell, G. Li . 
% In vivo contact kinematics and contact forces of the knee after total 
% knee  arthroplasty during dynamic weight-bearing activities. Journal of 
% Biomechanics 2008, 41(10):2159-68.
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
% Created by Rapha�l Dumas, Florent Moissent, Edouard Jouan
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
% April 2015
% Contact point constraints
%
% Modified by Raphael Dumas
% August 2018
% Renamed to multibody kinematics optimisation
%
% Modified by Raphael Dumas
% October 2019
% Patello-femoral mechanism personnalised with CT-scan geometry
% Renamed _GC5
% Possibility of informed segment length
%
% Modified by Raphael Dumas
% Januaray 2020
% Renamed _PS
%__________________________________________________________________________


function [Segment,Joint] = Multibody_Kinematics_Optimisation_PS(Segment,Joint,f)

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

% Mean segment geometry and markers coordinates
for i = 2:6 % From i = 2 (Foot) to i = 6 (Pelvis)
    
    % Segment length
    if ~(isfield(Segment,'L') && ~isempty(Segment(i).L)) % Informed segment length
        Segment(i).L = mean(sqrt(sum((Segment(i).Q(4:6,1,:) - ...
            Segment(i).Q(7:9,1,:)).^2)),3);  % Mean segment length
    end
    
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
rVs12 = [0; 0.3994 - 0.4260; 0]*(Segment(3).L/0.3994); % Scaled by shank length
% Shank length based on mid-maleolli (model virtual markers): 0.3994
% -0.4260 from knee joint centre to talus centre
% No longer need for rVs22, rVs32

% Ankle axis direction in shank segment
ns13 = [-0.105; -0.174; 0.979]/norm([-0.105; -0.174; 0.979]); % In Delp's thesis
% Ankle axis position in shank segment
% Witt regards to distal endpoints
% Ankle and subtalar axes superimposed (as well as Y-axes of foot and shank in neutral position)
rVs13 = rVs12;
% No longer need for rVs23, rVs33

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
% Contact points
% -------------------------------------------------------------------------
% Data from CAD of PS
% Component kinematics from Varadarajan et al. (2008)
% -------------------------------------------------------------------------

% Knee joint transformation matrix (from thigh to leg)
Joint(3).T = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(5).Q)),...
    Q2Tuv_array3(Segment(3).Q));
% Euler angles
Joint(3).Euler = R2mobileZXY_array3(Joint(3).T(1:3,1:3,:));
FE = permute(Joint(3).Euler(1,1,:),[3 1 2])/pi*180;% Extension(+)/Flexion(-)

% Polynomal
load  p_PS; % Femoral component and tibial insert are superimposed in neutral pose

% Virtual palpation in CT-scan/CAD model (tibia)
Tt_implant2bone = [ 0.9999    0.0112   -0.0093   -0.0027; ...
   -0.0111    0.9999    0.0105   -0.0043; ...
    0.0094   -0.0104    0.9999   -0.0027; ...
         0         0         0    1.0000];

% Virtual palpation in CT-scan/CAD model (femur)
Tf_implant2bone = [0.9999    0.0048   -0.0096   -0.0027; ...
   -0.0049    0.9999   -0.0135   -0.0043; ...
    0.0095    0.0135    0.9999   -0.0027; ...
         0         0         0    1.0000]; 

rVs43(:,1,:) = Mprod_array3(Tinv_array3(repmat(Tt_implant2bone, [1 1 n])), ... % Medial_contact_tibia (in m)
    [permute(polyval(pCM3x,FE),[2 3 1]); ...
    permute(polyval(pCM3y,FE),[2 3 1]); ...
    permute(polyval(pCM3z,FE),[2 3 1]); ...
    ones(1,1,n)]);
rVs43(4,:,:) = [];
rVs53(:,1,:) = Mprod_array3(Tinv_array3(repmat(Tt_implant2bone, [1 1 n])), ... % Lateral_contact_tibia (in m)
    [permute(polyval(pCL3x,FE),[2 3 1]); ... 
    permute(polyval(pCL3y,FE),[2 3 1]); ... 
    permute(polyval(pCL3z,FE),[2 3 1]); ...
    ones(1,1,n)]);
rVs53(4,:,:) = [];
% No longer need for rVs63, rVs73, rVs83

rVs15(:,1,:) = [0; - Segment(5).L; 0; 0] + ... Origin was at rD5
    Mprod_array3(Tinv_array3(repmat(Tf_implant2bone, [1 1 n])), ... % Medial_contact_femur (in m)
    [permute(polyval(pCM5x,FE),[2 3 1]); ...
    permute(polyval(pCM5y,FE),[2 3 1]); ...
    permute(polyval(pCM5z,FE),[2 3 1]); ...
    ones(1,1,n)]);
rVs15(4,:,:) = [];
rVs25(:,1,:) = [0; - Segment(5).L; 0; 0] + ... Origin was at rD5
    Mprod_array3(Tinv_array3(repmat(Tf_implant2bone, [1 1 n])), ... % Medial_contact_femur (in m)
    [permute(polyval(pCL5x,FE),[2 3 1]); ...
    permute(polyval(pCL5y,FE),[2 3 1]); ...
    permute(polyval(pCL5z,FE),[2 3 1]); ...
    ones(1,1,n)]);
rVs25(4,:,:) = [];
% No longer need for rVs35, rVs45, rVs55

% Coordinates of virtual markers in segment 3
% Expressed in (u3, rP3-rD3, w3)
nV43(:,:,:) = Mprod_array3(Minv_array3(repmat(Segment(3).B, [1,1,n])), rVs43); % V43: virtual marker 4 of segment 3
nV53(:,:,:) = Mprod_array3(Minv_array3(repmat(Segment(3).B, [1,1,n])), rVs53); % V53: virtual marker 5 of segment 3

% Coordinates of virtual markers in segment 5
% Expressed in (u5, rP5-rD5, w5)
nV15(:,:,:) = Mprod_array3(Minv_array3(repmat(Segment(5).B, [1,1,n])), rVs15); % V15: virtual marker 1 of segment 5
nV25(:,:,:) = Mprod_array3(Minv_array3(repmat(Segment(5).B, [1,1,n])), rVs25); % V25: virtual marker 2 of segment 5

% -------------------------------------------------------------------------
% Paralell patello-femoral mechanism: joint 4
% -------------------------------------------------------------------------
% Data from CT-scan and CAD of PS
% Axes of tibia, patella and femur assumed aligned in neutral posture
% -------------------------------------------------------------------------

% Hinge in femur (origin at distal endpoint)
rVs65 = [0; - Segment(5).L; 0; 0] + ... Origin was at rD5
    inv(Tf_implant2bone)*...
    [([-4.585; 4.718; 24.008*(-1)] + ...  % Medial sphere of femur throchlea (left knee)
    [-5.52; 3.136; -24.34*(-1)])... % Lateral sphere of femur throchlea (left knee)
    /2000;1]; % in m
rVs65(4,:) = [];

% Axis orientation
ns15 = ([-5.52; 3.136; -24.34*(-1)] - [-4.585; 4.718; 24.008*(-1)])/...
    norm([-5.52; 3.136; -24.34*(-1)] - [-4.585; 4.718; 24.008*(-1)]);

% Hinge in patella
rVs14 = [0; Segment(5).L; 0] + rVs65 - ...
    inv(Tf_implant2bone(1:3,1:3))*[42.1; -4.2; -2.72*(-1)]/1000; % Origin of patella (left knee)
% rVs24 = rD4:   Patellar tendon origin
ns14 = ns15; % Axes of patella and femur assumed aligned in neutral posture

% Patellar tendon insertion (origin at proximal endpoint)
rVs93 = inv(Tt_implant2bone)*...
    [[33.52;-79.87;-7.75*(-1)]/1000;1]; % Patellar tendon insertion (left knee)
rVs93(4,:) = [];

% Orientation of hinge axis
thetaP1 = acosd(dot(ns15,[1.0, 0.0, 0.0])); % Angle between n_axis_femur and X axis of femur SCS
thetaP2 = acosd(dot(ns15,[0.0, 1.0, 0.0])); % Angle between n_axis_femur and Y axis of femur SCS

% Patellar tendon length: distances associated with joint 4
Joint(4).d(1,1) = 55.127/1000;

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
Segment(5).nV(:,6) = inv(Segment(5).B)*rVs65; % V65: virtual marker 6 of segment 5
% Orientation (associated with segment 5)
Segment(5).nn(:,1) = inv(Segment(5).B)*ns15; % n15: virtual normal 1 of segment 5

% -------------------------------------------------------------------------
% Spherical hip joint: joint 5
% -------------------------------------------------------------------------
% Hip virtual marker mean coordinates (rV1 = rP6)
% Expressed in  in (u6, rP6-rD6, w6)
Segment(6).nV(:,1) = [0;-1;0] + inv(Segment(6).B)*...
    [0; 0; 0.0835]* ...% Hip joint centre on Z-axis
    (Segment(5).L/0.4040);  % Scaled on thigh length
% Thigh length based on mid-epicondyles (model virtual markers): 0.4040

% -------------------------------------------------------------------------
% Interpolation matrices
% -------------------------------------------------------------------------
% Segment 2
NV12 = [Segment(2).nV(1,1)*eye(3),...
    (1 + Segment(2).nV(2,1))*eye(3), ...
    - Segment(2).nV(2,1)*eye(3), ...
    Segment(2).nV(3,1)*eye(3)];
% NV22 = [Segment(2).nV(1,2)*eye(3),...
%     (1 + Segment(2).nV(2,2))*eye(3), ...
%     - Segment(2).nV(2,2)*eye(3), ...
%     Segment(2).nV(3,2)*eye(3)];

% Nn12 and Nn13 already computed


% Segment 3
NV13 = [Segment(3).nV(1,1)*eye(3),...
    (1 + Segment(3).nV(2,1))*eye(3), ...
    - Segment(3).nV(2,1)*eye(3), ...
    Segment(3).nV(3,1)*eye(3)];
Nn13 = [Segment(3).nn(1,1)*eye(3),...
    (Segment(3).nn(2,1))*eye(3), ...
    - Segment(3).nn(2,1)*eye(3), ...
    Segment(3).nn(3,1)*eye(3)];

NV43 = [nV43(1,1,:), zeros(1,2,n), ones(1,1,1) + nV43(2,1,:), zeros(1,2,n),  -nV43(2,1,:), zeros(1,2,n), nV43(3,1,:), zeros(1,2,n);...
    zeros(1,1,n), nV43(1,1,:), zeros(1,2,n), ones(1,1,1) + nV43(2,1,:), zeros(1,2,n),  -nV43(2,1,:), zeros(1,2,n), nV43(3,1,:), zeros(1,1,n); ...
    zeros(1,2,n), nV43(1,1,:), zeros(1,2,n), ones(1,1,1) + nV43(2,1,:), zeros(1,2,n),  -nV43(2,1,:), zeros(1,2,n), nV43(3,1,:)];
NV53 = [nV53(1,1,:), zeros(1,2,n), ones(1,1,1) + nV53(2,1,:), zeros(1,2,n),  -nV53(2,1,:), zeros(1,2,n), nV53(3,1,:), zeros(1,2,n);...
    zeros(1,1,n), nV53(1,1,:), zeros(1,2,n), ones(1,1,1) + nV53(2,1,:), zeros(1,2,n),  -nV53(2,1,:), zeros(1,2,n), nV53(3,1,:), zeros(1,1,n); ...
    zeros(1,2,n), nV53(1,1,:), zeros(1,2,n), ones(1,1,1) + nV53(2,1,:), zeros(1,2,n),  -nV53(2,1,:), zeros(1,2,n), nV53(3,1,:)];

NV93 = [Segment(3).nV(1,9)*eye(3),...
    (1 + Segment(3).nV(2,9))*eye(3), ...
    - Segment(3).nV(2,9)*eye(3), ...
    Segment(3).nV(3,9)*eye(3)];

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
NV15 = [nV15(1,1,:), zeros(1,2,n), ones(1,1,1) + nV15(2,1,:), zeros(1,2,n),  -nV15(2,1,:), zeros(1,2,n), nV15(3,1,:), zeros(1,2,n);...
    zeros(1,1,n), nV15(1,1,:), zeros(1,2,n), ones(1,1,1) + nV15(2,1,:), zeros(1,2,n),  -nV15(2,1,:), zeros(1,2,n), nV15(3,1,:), zeros(1,1,n); ...
    zeros(1,2,n), nV15(1,1,:), zeros(1,2,n), ones(1,1,1) + nV15(2,1,:), zeros(1,2,n),  -nV15(2,1,:), zeros(1,2,n), nV15(3,1,:)];
NV25 = [nV25(1,1,:), zeros(1,2,n), ones(1,1,1) + nV25(2,1,:), zeros(1,2,n),  -nV25(2,1,:), zeros(1,2,n), nV25(3,1,:), zeros(1,2,n);...
    zeros(1,1,n), nV25(1,1,:), zeros(1,2,n), ones(1,1,1) + nV25(2,1,:), zeros(1,2,n),  -nV25(2,1,:), zeros(1,2,n), nV25(3,1,:), zeros(1,1,n); ...
    zeros(1,2,n), nV25(1,1,:), zeros(1,2,n), ones(1,1,1) + nV25(2,1,:), zeros(1,2,n),  -nV25(2,1,:), zeros(1,2,n), nV25(3,1,:)];

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
    
    
    % Contact point trajectories
    % ---------------------------------------------------------------------
    % Position of virtual markers
    rV43 = Mprod_array3(NV43, Segment(3).Q);
    rV53 = Mprod_array3(NV53, Segment(3).Q);
    rV15 = Mprod_array3(NV15, Segment(5).Q);
    rV25 = Mprod_array3(NV25, Segment(5).Q);
    
    % Vector of kinematic constraints
    phikK = [rV15 - rV43; ... % Medial contact
    rV25(1:2,1,:) - rV53(1:2,1,:)]; % Lateral contact

   % Jacobian of kinematic constraints
    KkK = zeros(5,5*12,n); % Initialisation
    KkK(1:3,13:24,:) = - NV43;
    KkK(1:3,37:48,:) = NV15;
    KkK(4:5,13:24,:) = - NV53(1:2,:,:);
    KkK(4:5,37:48,:) = NV25(1:2,:,:);
    Joint(3).Kk = KkK; 
    
    % Partial derivative of Jacobian * Lagrange multipliers
    dKlambdakKdQ = zeros(5*12,5*12,n); % Initialisation
    
    Segment(3).rV43 = Mprod_array3(NV43, Segment(3).Q); % Stored in Segment
    Segment(3).rV53 = Mprod_array3(NV53, Segment(3).Q); % Stored in Segment
    Segment(5).rV15 = Mprod_array3(NV15, Segment(5).Q); % Stored in Segment
    Segment(5).rV25 = Mprod_array3(NV25, Segment(5).Q); % Stored in Segment

    
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
% Inconsistent velocities
% -------------------------------------------------------------------------
% Q
Q = [Segment(2).Q; ...   % Foot
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

% Tibio-femoral
% X3
NX3 = [eye(3),zeros(3,9)]; % X3 = u3
X3 = Segment(3).Q(1:3,1,:); % u3 = X3 of SCS
% Y3
inB3 = inv(Segment(3).B);
NY3 = [inB3(1,2)*eye(3),inB3(2,2)*eye(3),-inB3(2,2)*eye(3),inB3(3,2)*eye(3)];
Y3 = Mprod_array3(repmat(NY3,[1,1,n]),Segment(3).Q); % Y3 of SCS
% Z3
NZ3 = [inB3(1,3)*eye(3),inB3(2,3)*eye(3),-inB3(2,3)*eye(3),inB3(3,3)*eye(3)];
Z3 = Mprod_array3(repmat(NZ3,[1,1,n]),Segment(3).Q); % Z3 of SCS
%
Joint(3).Kk(1,13:24,:) = -Mprod_array3(permute(X3,[2,1,3]), NV43) + ...
    Mprod_array3(permute(rV15 - rV43,[2,1,3]),repmat(NX3,[1,1,n]));
Joint(3).Kk(1,37:48,:) = Mprod_array3(permute(X3,[2,1,3]), NV15);
Joint(3).Kk(2,13:24,:) = -Mprod_array3(permute(Y3,[2,1,3]), NV43) + ...
    Mprod_array3(permute(rV15 - rV43,[2,1,3]),repmat(NY3,[1,1,n]));
Joint(3).Kk(2,37:48,:) = Mprod_array3(permute(Y3,[2,1,3]), NV15);
Joint(3).Kk(3,13:24,:) = -Mprod_array3(permute(Z3,[2,1,3]), NV43) + ...
    Mprod_array3(permute(rV15 - rV43,[2,1,3]),repmat(NZ3,[1,1,n]));
Joint(3).Kk(3,37:48,:) = Mprod_array3(permute(Z3,[2,1,3]), NV15);
Joint(3).Kk(4,13:24,:) = -Mprod_array3(permute(X3,[2,1,3]), NV53) + ...
    Mprod_array3(permute(rV25 - rV53,[2,1,3]),repmat(NX3,[1,1,n]));
Joint(3).Kk(4,37:48,:) = Mprod_array3(permute(X3,[2,1,3]), NV25);
Joint(3).Kk(5,13:24,:) = -Mprod_array3(permute(Y3,[2,1,3]), NV53) + ...
    Mprod_array3(permute(rV25 - rV53,[2,1,3]),repmat(NY3,[1,1,n]));
Joint(3).Kk(5,37:48,:) = Mprod_array3(permute(Y3,[2,1,3]), NV25);

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
