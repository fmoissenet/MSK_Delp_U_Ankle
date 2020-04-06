% FUNCTION
% Multibody_Kinematics_Optimisation.m
%__________________________________________________________________________
%
% PURPOSE
% Multibody optimisation of the postions of foot, shank, patella, thigh
% and computation of their consistent velocities and accelerations with
% universal ankle, tibio-femoral contact points, paralell patello-femoral 
% mechanism and spherial hip joint
%
% SYNOPSIS
% [Segment,Joint] = Multibody_Kinematics_Optimisation(Segment,Joint,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Joint (cf. data structure in user guide)
% Model to pass number of frames, acquisition and cutt off frequencies
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
% System Dynamics 2012, 28:125–141.
% S. Duprey, L. Cheze, R. Dumas. Influence of joint constraints on lower 
% limb kinematics estimation from skin markers using global optimization.
% Journal of Biomechancis 2010, 43(14):2858-62.
% A. Zeighami, R. Aissaoui, R. Dumas. Knee medial and lateral contact 
% forces in a musculoskeletal model with subject-specific contact point 
% trajectories. Journal of Biomechanics 2018, 69: 138–145.
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Vnop_array3.m
% Mprod_array3.m
% Derive_array3.m
% Vfilt_array3.m
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
% September 2019
% Patella hinge based on SIMM data
%
% Modified by Raphael Dumas
% March 2020
% Generic vs. Informed structures (Segment, Joint, Model)
% n,f,fc as Model.Informed fields
%__________________________________________________________________________


function [Segment,Joint] = Multibody_Kinematics_Optimisation(Segment,Joint,Model)

% Number of frames
n = Model.Informed.n;
% Acquisition frequency
f = Model.Informed.f;
% Cut off frequency for filtering
fc = Model.Informed.fc;


%% -------------------------------------------------------------------------
% Segment and Joint parameters
% -------------------------------------------------------------------------

% Initialisation
Segment(1).Informed.L = NaN; % No value for segment 1 (Forceplate)
Segment(1).Informed.alpha = NaN; % No value for segment 1 (Forceplate)
Segment(1).Informed.beta = NaN; % No value for segment 1 (Forceplate)
Segment(1).Informed.gamma = NaN; % No value for segment 1 (Forceplate)

% Mean segment geometry and markers coordinates
for i = 2:6 % From i = 2 (Foot) to i = 6 (Pelvis)
    % Mean coordinates of markers in (u, rP-rD, w)
    for j = 1:size(Segment(i).Informed.rM,2)
        % Projection in a non orthonormal coordinate system
        Segment(i).Informed.nM(:,j) = mean(Vnop_array3(...
            Segment(i).Informed.rM(:,j,:) - Segment(i).Informed.Q(4:6,1,:),...
            Segment(i).Informed.Q(1:3,1,:),...
            Segment(i).Informed.Q(4:6,1,:) - Segment(i).Informed.Q(7:9,1,:),...
            Segment(i).Informed.Q(10:12,1,:)),3);
    end
end

% Initialisation
Joint(1).Informed.Kk = [];


%% -------------------------------------------------------------------------
% Universal ankle joint
% -------------------------------------------------------------------------

% Coordinates of virtual markers in Segment 2
% Expressed in (u2, rP2-rD2, w2)
% Subtalar axis position in foot segment
% Scaled by a homothety of centre at proximal endpoint
Segment(2).Informed.nV(:,1) = inv(Segment(2).Informed.B) * ...
    Segment(2).Generic.rVs(:,1)*Segment(2).Informed.Scale; % Virtual marker 1 of segment 2
% Subtalar axis direction in foot segment
Segment(2).Informed.nn(:,1) = inv(Segment(2).Informed.B) * ...
    Segment(2).Generic.ns(:,1); % Virtual direction 1 of segment 2

% Coordinates of virtual markers in Segment 3
% Expressed in (u3, rP3-rD3, w3)
% Ankle axis position in shank segment
% Scaled by a homothety of centre at distal endpoint
Segment(3).Informed.nV(:,1) = [0; -1; 0] + inv(Segment(3).Informed.B) * ...
        (Segment(3).Generic.rVs(:,1) + [0; Segment(3).Generic.L; 0])*...
        Segment(3).Informed.Scale;
% Ankle axis direction in shank segment
Segment(3).Informed.nn(:,1) = inv(Segment(3).Informed.B) * ...
    Segment(3).Generic.ns(:,1); % Virtual direction 1 of segment 3

% Angle between universal joint axes
Nn12 = [Segment(2).Informed.nn(1,1)*eye(3),...
    (Segment(2).Informed.nn(2,1))*eye(3), ...
    - Segment(2).Informed.nn(2,1)*eye(3), ...
    Segment(2).Informed.nn(3,1)*eye(3)]; % Interpolation matrix
Nn13 = [Segment(3).Informed.nn(1,1)*eye(3),...
    (Segment(3).Informed.nn(2,1))*eye(3), ...
    - Segment(3).Informed.nn(2,1)*eye(3), ...
    Segment(3).Informed.nn(3,1)*eye(3)]; % Interpolation matrix

% Interpolation matrices for Segment 2
% Nn12 already computed
NV12 = [Segment(2).Informed.nV(1,1)*eye(3),...
    (1 + Segment(2).Informed.nV(2,1))*eye(3), ...
    - Segment(2).Informed.nV(2,1)*eye(3), ...
    Segment(2).Informed.nV(3,1)*eye(3)];

% Interpolation matrices for Segment 3
% Nn13 already computed
NV13 = [Segment(3).Informed.nV(1,1)*eye(3),...
    (1 + Segment(3).Informed.nV(2,1))*eye(3), ...
    - Segment(3).Informed.nV(2,1)*eye(3), ...
    Segment(3).Informed.nV(3,1)*eye(3)];


%% -------------------------------------------------------------------------
% Tibiofemoral contact points trajectories
% -------------------------------------------------------------------------

% Coordinates of virtual markers in Segment 3
% Expressed in (u3, rP3-rD3, w3)
% Medial contact in tibia
% Scaled by a homothety of centre at proximal endpoint
nV43(:,:,:) = Mprod_array3(Minv_array3(repmat(Segment(3).Informed.B,[1,1,n])), ...
    Segment(3).Generic.rVs43*Segment(3).Informed.Scale); % Virtual marker 4 of segment 3 (depending on time)
% Lateral contact in tibia
nV53(:,:,:) = Mprod_array3(Minv_array3(repmat(Segment(3).Informed.B,[1,1,n])), ...
    Segment(3).Generic.rVs53*Segment(3).Informed.Scale); % Virtual marker 5 of segment 3 (depending on time)
% No need for other virtual markers (rVs63, rVs73, and rVs83)

% Coordinates of virtual markers in Segment 5
% Expressed in (u5, rP5-rD5, w5)
% Medial contact in femur
% Scaled by a homothety of centre at distal endpoint
nV15(:,:,:) = repmat([0; -1; 0],[1,1,n]) + ...
    Mprod_array3(Minv_array3(repmat(Segment(5).Informed.B,[1,1,n])), ...
    (Segment(5).Generic.rVs15 + repmat([0; Segment(5).Generic.L; 0],[1,1,n])) * ...
    Segment(5).Informed.Scale); % Virtual marker 1 of segment 5 (depending on time)
% Lateral contact in femur
nV25(:,:,:) = repmat([0; -1; 0],[1,1,n]) + ...
    Mprod_array3(Minv_array3(repmat(Segment(5).Informed.B,[1,1,n])), ...
    (Segment(5).Generic.rVs25 + repmat([0; Segment(5).Generic.L; 0],[1,1,n])) * ...
    Segment(5).Informed.Scale); % Virtual marker 2 of segment 5 (depending on time)
% No need for other virtual markers (rVs35, rVs45, and rVs55)

% Interpolation matrices for Segment 3
NV43 = [nV43(1,1,:), zeros(1,2,n), ones(1,1,n) + nV43(2,1,:), zeros(1,2,n), -nV43(2,1,:), zeros(1,2,n), nV43(3,1,:), zeros(1,2,n);...
    zeros(1,1,n), nV43(1,1,:), zeros(1,2,n), ones(1,1,n) + nV43(2,1,:), zeros(1,2,n), -nV43(2,1,:), zeros(1,2,n), nV43(3,1,:), zeros(1,1,n); ...
    zeros(1,2,n), nV43(1,1,:), zeros(1,2,n), ones(1,1,n) + nV43(2,1,:), zeros(1,2,n), -nV43(2,1,:), zeros(1,2,n), nV43(3,1,:)];
NV53 = [nV53(1,1,:), zeros(1,2,n), ones(1,1,n) + nV53(2,1,:), zeros(1,2,n), -nV53(2,1,:), zeros(1,2,n), nV53(3,1,:), zeros(1,2,n);...
    zeros(1,1,n), nV53(1,1,:), zeros(1,2,n), ones(1,1,n) + nV53(2,1,:), zeros(1,2,n), -nV53(2,1,:), zeros(1,2,n), nV53(3,1,:), zeros(1,1,n); ...
    zeros(1,2,n), nV53(1,1,:), zeros(1,2,n), ones(1,1,n) + nV53(2,1,:), zeros(1,2,n), -nV53(2,1,:), zeros(1,2,n), nV53(3,1,:)];

% Interpolation matrices for Segment 5
NV15 = [nV15(1,1,:), zeros(1,2,n), ones(1,1,n) + nV15(2,1,:), zeros(1,2,n), -nV15(2,1,:), zeros(1,2,n), nV15(3,1,:), zeros(1,2,n);...
    zeros(1,1,n), nV15(1,1,:), zeros(1,2,n), ones(1,1,n) + nV15(2,1,:), zeros(1,2,n), -nV15(2,1,:), zeros(1,2,n), nV15(3,1,:), zeros(1,1,n); ...
    zeros(1,2,n), nV15(1,1,:), zeros(1,2,n), ones(1,1,n) + nV15(2,1,:), zeros(1,2,n), -nV15(2,1,:), zeros(1,2,n), nV15(3,1,:)];
NV25 = [nV25(1,1,:), zeros(1,2,n), ones(1,1,n) + nV25(2,1,:), zeros(1,2,n), -nV25(2,1,:), zeros(1,2,n), nV25(3,1,:), zeros(1,2,n);...
    zeros(1,1,n), nV25(1,1,:), zeros(1,2,n), ones(1,1,n) + nV25(2,1,:), zeros(1,2,n), -nV25(2,1,:), zeros(1,2,n), nV25(3,1,:), zeros(1,1,n); ...
    zeros(1,2,n), nV25(1,1,:), zeros(1,2,n), ones(1,1,n) + nV25(2,1,:), zeros(1,2,n), -nV25(2,1,:), zeros(1,2,n), nV25(3,1,:)];


%% -------------------------------------------------------------------------
% ¨Parallel patello-femoral joint
% -------------------------------------------------------------------------

% Coordinates of virtual markers in Segment 3
% Expressed in (u3, rP3-rD3, w3)
% Patellar tendon insertion
% Scaled by a homothety of centre at proximal endpoint
Segment(3).Informed.nV(:,9) = inv(Segment(3).Informed.B) *...
    Segment(3).Generic.rVs(:,9)*Segment(3).Informed.Scale; % V93: virtual marker 9 of segment 3

% Coordinates of virtual markers in Segment 4
% Expressed in (u4, rP4-rD4, w4)
% Position of hinge axis
Segment(4).Informed.nV(:,1) = inv(Segment(4).Informed.B)* ...
    Segment(4).Generic.rVs(:,1)*Segment(4).Informed.Scale; % V14: virtual marker 1 of segment 4
% Orientation of hinge axis
Segment(4).Informed.nn(:,1) = inv(Segment(4).Informed.B) * ...
    Segment(4).Generic.ns(:,1); % n14: virtual direction 1 of segment 4

% Coordinates of virtual markers in Segment 5
% Expressed in (u5, rP5-rD5, w5)
% Position of hinge axis
% Scaled by a homothety of centre at distal endpoint
Segment(5).Informed.nV(:,6) = [0;-1;0] + ...
    inv(Segment(5).Informed.B) * ...
    (Segment(5).Generic.rVs(:,6) + [0; Segment(5).Generic.L; 0]) * ...
    Segment(5).Informed.Scale; % V65: virtual marker 6 of segment 5
% Orientation of hinge axis
Segment(5).Informed.nn(:,1) = inv(Segment(5).Informed.B) * ...
    Segment(5).Generic.ns(:,1); % n15: virtual direction 1 of segment 5

% Interpolation matrix for Segment 3
NV93 = [Segment(3).Informed.nV(1,9)*eye(3),...
    (1 + Segment(3).Informed.nV(2,9))*eye(3), ...
    - Segment(3).Informed.nV(2,9)*eye(3), ...
    Segment(3).Informed.nV(3,9)*eye(3)];

% Interpolation matrices for Segment 3
NV14 = [Segment(4).Informed.nV(1,1)*eye(3),...
    (1 + Segment(4).Informed.nV(2,1))*eye(3), ...
    - Segment(4).Informed.nV(2,1)*eye(3), ...
    Segment(4).Informed.nV(3,1)*eye(3)];
Nn14 = [Segment(4).Informed.nn(1,1)*eye(3),...
    (Segment(4).Informed.nn(2,1))*eye(3), ...
    - Segment(4).Informed.nn(2,1)*eye(3), ...
    Segment(4).Informed.nn(3,1)*eye(3)];

% Interpolation matrices for Segment 5
NV65 = [Segment(5).Informed.nV(1,6)*eye(3),...
    (1 + Segment(5).Informed.nV(2,6))*eye(3), ...
    - Segment(5).Informed.nV(2,6)*eye(3), ...
    Segment(5).Informed.nV(3,6)*eye(3)];


%% -------------------------------------------------------------------------
% Spherical hip joint
% -------------------------------------------------------------------------

% Coordinates of virtual markers in Segment 6
% Expressed in  in (u6, rP6-rD6, w6)
% Scaled by a homothety of centre at proximal endpoint
Segment(6).Informed.nV(:,1) =  inv(Segment(6).Informed.B) * ...
    Segment(6).Generic.rVs(:,1)*...
    Segment(6).Informed.Scale;

% Interpolation matrix for Segment 6
NV16 = [Segment(6).Informed.nV(1,1)*eye(3),...
    (1 + Segment(6).Informed.nV(2,1))*eye(3), ...
    - Segment(6).Informed.nV(2,1)*eye(3), ...
    Segment(6).Informed.nV(3,1)*eye(3)];


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
    rV12 = Mprod_array3(repmat(NV12,[1,1,n]),Segment(2).Informed.Q); % Position of subtalar axis
    rV13 = Mprod_array3(repmat(NV13,[1,1,n]),Segment(3).Informed.Q); % Position of ankle axis
    % Direction of virtual orientation
    n12 = Mprod_array3(repmat(Nn12,[1,1,n]),Segment(2).Informed.Q); % Orientation of subtalar axis
    n13 = Mprod_array3(repmat(Nn13,[1,1,n]),Segment(3).Informed.Q); % Orientation of ankle axis
    % Angle
    thetaA = mean(acosd(dot(n13,n12)),3); % Mean angle between subtalar and ankle axes 
        
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
    Joint(2).Informed.Kk = KkA;
    
    % Partial derivative of Jacobian * Lagrange multipliers
    dKlambdakAdQ = zeros(5*12,5*12,n); % Initialisation
    dKlambdakAdQ(1:12,13:24,:) = Mprod_array3(lambdakA(4,1,:),repmat(Nn12'*Nn13,[1,1,n]));
    dKlambdakAdQ(13:24,1:12,:) = Mprod_array3(lambdakA(4,1,:),repmat(Nn13'*Nn12,[1,1,n]));
    
    % Segment structure
    Segment(2).Informed.n12 = n12; % To be potentially plot
    Segment(3).Informed.n13 = n13; % To be potentially plot

    
    % Tibiofemoral contact point trajectories
    % ---------------------------------------------------------------------
    
    % Position of virtual markers
    rV43 = Mprod_array3(NV43,Segment(3).Informed.Q); 
    rV53 = Mprod_array3(NV53,Segment(3).Informed.Q); 
    rV15 = Mprod_array3(NV15,Segment(5).Informed.Q); 
    rV25 = Mprod_array3(NV25,Segment(5).Informed.Q);
    
    % Vector of kinematic constraints
    phikK = [rV15 - rV43; ... % Medial contact
    rV25(1:2,1,:) - rV53(1:2,1,:)]; % Lateral contact

   % Jacobian of kinematic constraints
    KkK = zeros(5,5*12,n); % Initialisation
    KkK(1:3,13:24,:) = - NV43;
    KkK(1:3,37:48,:) = NV15;
    KkK(4:5,13:24,:) = - NV53(1:2,:,:);
    KkK(4:5,37:48,:) = NV25(1:2,:,:);
    
    % Joint structure
    Joint(3).Informed.Kk = KkK; 
    
    % Partial derivative of Jacobian * Lagrange multipliers
    dKlambdakKdQ = zeros(5*12,5*12,n); % Initialisation
    
    % Segment structure
    Segment(3).Informed.rV43 = rV43; % To be potentially plot
    Segment(3).Informed.rV53 = rV53; % To be potentially plot
    Segment(5).Informed.rV15 = rV15; % To be potentially plot
    Segment(5).Informed.rV25 = rV25; % To be potentially plot
     
    % Parallel patello-femoral mechanism
    % ---------------------------------------------------------------------
    % Position of virtual markers
    rV93 = Mprod_array3(repmat(NV93,[1,1,n]),Segment(3).Informed.Q); % Patellar tendon insertion
    rV14 = Mprod_array3(repmat(NV14,[1,1,n]),Segment(4).Informed.Q); % Position of hinge axis
    rV65 = Mprod_array3(repmat(NV65,[1,1,n]),Segment(5).Informed.Q); % Position of hinge axis
    % Orientation of virtual direction
    n14 = Mprod_array3(repmat(Nn14,[1,1,n]),Segment(4).Informed.Q);  % Orientation of hinge axis
    % Angles
    thetaP1 = acosd(dot(Segment(5).Informed.nn(:,1),[1; 0; 0])); % Between hinge axis and X-axis of Thigh
    thetaP2 = acosd(dot(Segment(5).Informed.nn(:,1),[0; 1; 0])); % Between hinge axis and Y-axis of Thigh 
    % Patellar tendon length
    Joint(4).Informed.d(1,1) = norm([0; Segment(5).Informed.L; 0] ... % From rD to rP in Thigh
        + Segment(5).Generic.rVs(:,6)*Segment(5).Informed.Scale ... % From rP to position of hinge axis in Thigh
        - Segment(4).Generic.rVs(:,1)*Segment(4).Informed.Scale ...  % From position of hinge axis to rP in Patella
        - [0; Segment(4).Generic.L; 0]*Segment(4).Informed.Scale ... % From rP to rD in Patella
        - Segment(3).Generic.rVs(:,9)*Segment(3).Informed.Scale); % From patellar tendon insertion to rP in Shank (= rD in Thigh)

    % Vector of kinematic constraints
    phikP = [rV65 - rV14; ...
        dot(Segment(5).Informed.Q(1:3,1,:),n14)- cosd(thetaP1); ...
        dot((Segment(5).Informed.Q(4:6,1,:)- Segment(5).Informed.Q(7:9,1,:)), n14) - ...
        Segment(5).Informed.L*cosd(thetaP2); ...
        dot((Segment(4).Informed.Q(7:9,1,:) - rV93),(Segment(4).Informed.Q(7:9,1,:) - rV93)) - ...
        (repmat(Joint(4).Informed.d(1,1)^2,[1,1,n]))];
    
    
    % Jacobian of kinematic constraints
    KkP = zeros(6,5*12,n); % Initialisation
    KkP(1:3,25:36,:) = repmat(-NV14,[1,1,n]);
    KkP(1:3,37:48,:) = repmat(NV65,[1,1,n]);
    KkP(4,25:36,:) = Mprod_array3(permute(Segment(5).Informed.Q(1:3,1,:),[2,1,3]),repmat(Nn14,[1,1,n]));
    KkP(4,37:48,:) = [permute(n14,[2,1,3]),zeros(1,3,n),zeros(1,3,n),zeros(1,3,n)];
    KkP(5,25:36,:) = Mprod_array3(permute(Segment(5).Informed.Q(4:6,1,:) - Segment(5).Informed.Q(7:9,1,:),[2,1,3]),repmat(Nn14,[1,1,n]));
    KkP(5,37:48,:) = [zeros(1,3,n),permute(n14,[2,1,3]),permute(-n14,[2,1,3]),zeros(1,3,n)];
    KkP(6,13:24,:) = - Mprod_array3(permute(Segment(4).Informed.Q(7:9,1,:) - rV93,[2,1,3]),repmat(2*NV93,[1,1,n]));
    KkP(6,31:33,:) = Mprod_array3(permute(Segment(4).Informed.Q(7:9,1,:) - rV93,[2,1,3]),repmat(2*eye(3),[1,1,n]));
    % with transpose = permute( ,[2,1,3])
    
    % Joint structure
    Joint(4).Informed.Kk = KkP;
    
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
    
    % Segment structure
    Segment(4).Informed.n14 = n14; % To be potentially plot
    
    % Spherical hip joint
    % ---------------------------------------------------------------------
    
    % Position of virtual markers
    rV16 = Mprod_array3(repmat(NV16,[1,1,n]),Segment(6).Informed.Q); % HJC_pelvis
    
    % Vector of kinematic constraints
    phikH = rV16 - Segment(5).Informed.Q(4:6,1,:);
    
    % Jacobian of kinematic constraints
    KkH = zeros(3,5*12,n); % Initialisation
    KkH(1:3,40:42,:) = repmat(-eye(3),[1,1,n]);
    KkH(1:3,49:60,:) = repmat(NV16,[1,1,n]);
    
    % Joint structure
    Joint(5).Informed.Kk = KkH;
    
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
        ui = Segment(i).Informed.Q(1:3,1,:);
        vi = Segment(i).Informed.Q(4:6,1,:) - Segment(i).Informed.Q(7:9,1,:);
        wi = Segment(i).Informed.Q(10:12,1,:);
        phiri = [dot(ui,ui) - ones(1,1,n);...
            dot(ui,vi) - repmat(Segment(i).Informed.L*cosd(Segment(i).Informed.gamma),[1,1,n]); ...
            dot(ui,wi) - repmat(cosd(Segment(i).Informed.beta),[1,1,n]); ...
            dot(vi,vi) - repmat(Segment(i).Informed.L^2,[1,1,n]);
            dot(vi,wi) - repmat(Segment(i).Informed.L*cosd(Segment(i).Informed.alpha),[1,1,n]);
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
        Segment(i).Informed.Kr = Kri;
        
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
        Kmi = zeros(size(Segment(i).Informed.rM,2)*3,5*12,n); % Initialisation
        phimi = []; % Initialisation
        
        for j = 1:size(Segment(i).Informed.rM,2)
            % Interpolation matrix
            NMij = [Segment(i).Informed.nM(1,j)*eye(3),...
                (1 + Segment(i).Informed.nM(2,j))*eye(3), ...
                - Segment(i).Informed.nM(2,j)*eye(3), ...
                Segment(i).Informed.nM(3,j)*eye(3)];
            % Vector of driving constraints
            phimi((j-1)*3+1:(j-1)*3+3,1,:) = Segment(i).Informed.rM(:,j,:) ...
                - Mprod_array3(repmat(NMij,[1,1,n]),Segment(i).Informed.Q);
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
    Segment(2).Informed.Q = Segment(2).Informed.Q + dX(1:12,1,:);     % Foot
    Segment(3).Informed.Q = Segment(3).Informed.Q + dX(13:24,1,:);    % Shank
    Segment(4).Informed.Q = Segment(4).Informed.Q + dX(25:36,1,:);    % Patella
    Segment(5).Informed.Q = Segment(5).Informed.Q + dX(37:48,1,:);    % Thigh
    Segment(6).Informed.Q = Segment(6).Informed.Q + dX(49:60,1,:);    % Pelvis
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
Q = [Segment(2).Informed.Q; ...   % Foot
    Segment(3).Informed.Q; ...   % Shank
    Segment(4).Informed.Q; ...   % Patella
    Segment(5).Informed.Q; ...   % Thigh
    Segment(6).Informed.Q];      % Pelvis

% dQdt with filtering
dQdt = Vfilt_array3(Derive_array3(Q,1/f),f,fc);
% Extraction from dQdt
Segment(2).Informed.dQdt = dQdt(1:12,1,:); % Foot
Segment(3).Informed.dQdt = dQdt(13:24,1,:); % Shank
Segment(4).Informed.dQdt = dQdt(25:36,1,:); % Patella
Segment(5).Informed.dQdt = dQdt(37:48,1,:); % Thigh
Segment(6).Informed.dQdt = dQdt(49:60,1,:); % Pelvis

% d2Qdt2 with filtering
d2Qdt2 = Vfilt_array3(Derive_array3(dQdt,1/f),f,fc);
% Extraction from d2Qdt2
Segment(2).Informed.d2Qdt2 = d2Qdt2(1:12,1,:); % Foot
Segment(3).Informed.d2Qdt2 = d2Qdt2(13:24,1,:); % Shank
Segment(4).Informed.d2Qdt2 = d2Qdt2(25:36,1,:); % Patella
Segment(5).Informed.d2Qdt2 = d2Qdt2(37:48,1,:); % Thigh
Segment(6).Informed.d2Qdt2 = d2Qdt2(49:60,1,:); % Pelvis


%% -------------------------------------------------------------------------
% Modified Jacobian of the kinematic constraint for spherical joints
% Constraints projected on the axes of the distal segment
% -------------------------------------------------------------------------

% Ankle
% X2
NX2 = [eye(3),zeros(3,9)]; % X2 = u2
X2 = Segment(2).Informed.Q(1:3,1,:); % u2 = X2 of SCS
% Y2
inB2 = inv(Segment(2).Informed.B);
NY2 = [inB2(1,2)*eye(3),inB2(2,2)*eye(3),-inB2(2,2)*eye(3),inB2(3,2)*eye(3)];
Y2 = Mprod_array3(repmat(NY2,[1,1,n]),Segment(2).Informed.Q); % Y2 of SCS
% Z2
NZ2 = [inB2(1,3)*eye(3),inB2(2,3)*eye(3),-inB2(2,3)*eye(3),inB2(3,3)*eye(3)];
Z2 = Mprod_array3(repmat(NZ2,[1,1,n]),Segment(2).Informed.Q); % Z2 of SCS
%
Joint(2).Informed.Kk(1,1:12,:) = -Mprod_array3(permute(X2,[2,1,3]),repmat(NV12,[1,1,n])) + ...
    Mprod_array3(permute(rV13 - rV12,[2,1,3]),repmat(NX2,[1,1,n]));
Joint(2).Informed.Kk(1,13:24,:) = Mprod_array3(permute(X2,[2,1,3]),repmat(NV13,[1,1,n]));
Joint(2).Informed.Kk(2,1:12,:) = -Mprod_array3(permute(Y2,[2,1,3]),repmat(NV12,[1,1,n])) + ...
    Mprod_array3(permute(rV13 - rV12,[2,1,3]),repmat(NY2,[1,1,n]));
Joint(2).Informed.Kk(2,13:24,:) = Mprod_array3(permute(Y2,[2,1,3]),repmat(NV13,[1,1,n]));
Joint(2).Informed.Kk(3,1:12,:) = -Mprod_array3(permute(Z2,[2,1,3]),repmat(NV12,[1,1,n])) + ...
    Mprod_array3(permute(rV13 - rV12,[2,1,3]),repmat(NZ2,[1,1,n]));
Joint(2).Informed.Kk(3,13:24,:) = Mprod_array3(permute(Z2,[2,1,3]),repmat(NV13,[1,1,n]));

% Tibio-femoral
% X3
NX3 = [eye(3),zeros(3,9)]; % X3 = u3
X3 = Segment(3).Informed.Q(1:3,1,:); % u3 = X3 of SCS
% Y3
inB3 = inv(Segment(3).Informed.B);
NY3 = [inB3(1,2)*eye(3),inB3(2,2)*eye(3),-inB3(2,2)*eye(3),inB3(3,2)*eye(3)];
Y3 = Mprod_array3(repmat(NY3,[1,1,n]),Segment(3).Informed.Q); % Y3 of SCS
% Z3
NZ3 = [inB3(1,3)*eye(3),inB3(2,3)*eye(3),-inB3(2,3)*eye(3),inB3(3,3)*eye(3)];
Z3 = Mprod_array3(repmat(NZ3,[1,1,n]),Segment(3).Informed.Q); % Z3 of SCS
%
Joint(3).Informed.Kk(1,13:24,:) = -Mprod_array3(permute(X3,[2,1,3]),NV43) + ...
    Mprod_array3(permute(rV15 - rV43,[2,1,3]),repmat(NX3,[1,1,n]));
Joint(3).Informed.Kk(1,37:48,:) = Mprod_array3(permute(X3,[2,1,3]),NV15);
Joint(3).Informed.Kk(2,13:24,:) = -Mprod_array3(permute(Y3,[2,1,3]),NV43) + ...
    Mprod_array3(permute(rV15 - rV43,[2,1,3]),repmat(NY3,[1,1,n]));
Joint(3).Informed.Kk(2,37:48,:) = Mprod_array3(permute(Y3,[2,1,3]),NV15);
Joint(3).Informed.Kk(3,13:24,:) = -Mprod_array3(permute(Z3,[2,1,3]),NV43) + ...
    Mprod_array3(permute(rV15 - rV43,[2,1,3]),repmat(NZ3,[1,1,n]));
Joint(3).Informed.Kk(3,37:48,:) = Mprod_array3(permute(Z3,[2,1,3]),NV15);
Joint(3).Informed.Kk(4,13:24,:) = -Mprod_array3(permute(X3,[2,1,3]),NV53) + ...
    Mprod_array3(permute(rV25 - rV53,[2,1,3]),repmat(NX3,[1,1,n]));
Joint(3).Informed.Kk(4,37:48,:) = Mprod_array3(permute(X3,[2,1,3]),NV25);
Joint(3).Informed.Kk(5,13:24,:) = -Mprod_array3(permute(Y3,[2,1,3]),NV53) + ...
    Mprod_array3(permute(rV25 - rV53,[2,1,3]),repmat(NY3,[1,1,n]));
Joint(3).Informed.Kk(5,37:48,:) = Mprod_array3(permute(Y3,[2,1,3]),NV25);

% Patello-femoral
% X4
NX4 = [eye(3),zeros(3,9)]; % X4 = u4
X4 = Segment(4).Informed.Q(1:3,1,:); % u4 = X4 of SCS
% Y4
inB4 = inv(Segment(4).Informed.B);
NY4 = [inB4(1,2)*eye(3),inB4(2,2)*eye(3),-inB4(2,2)*eye(3),inB4(3,2)*eye(3)];
Y4 = Mprod_array3(repmat(NY4,[1,1,n]),Segment(4).Informed.Q); % Y4 of SCS
% Z4
NZ4 = [inB4(1,3)*eye(3),inB4(2,3)*eye(3),-inB4(2,3)*eye(3),inB4(3,3)*eye(3)];
Z4 = Mprod_array3(repmat(NZ4,[1,1,n]),Segment(4).Informed.Q); % Z4 of SCS
%
Joint(4).Informed.Kk(1,25:36,:) = -Mprod_array3(permute(X4,[2,1,3]),repmat(NV14,[1,1,n])) + ...
    Mprod_array3(permute(rV65 - rV14,[2,1,3]),repmat(NX4,[1,1,n]));
Joint(4).Informed.Kk(1,37:48,:) = Mprod_array3(permute(X4,[2,1,3]),repmat(NV65,[1,1,n]));
Joint(4).Informed.Kk(2,25:36,:) = -Mprod_array3(permute(Y4,[2,1,3]),repmat(NV14,[1,1,n])) + ...
    Mprod_array3(permute(rV65 - rV14,[2,1,3]),repmat(NY4,[1,1,n]));
Joint(4).Informed.Kk(2,37:48,:) = Mprod_array3(permute(Y4,[2,1,3]),repmat(NV65,[1,1,n]));
Joint(4).Informed.Kk(3,25:36,:) = -Mprod_array3(permute(Z4,[2,1,3]),repmat(NV14,[1,1,n])) + ...
    Mprod_array3(permute(rV65 - rV14,[2,1,3]),repmat(NZ4,[1,1,n]));
Joint(4).Informed.Kk(3,37:48,:) = Mprod_array3(permute(Z4,[2,1,3]),repmat(NV65,[1,1,n]));

% Hip
% X5
NX5 = [eye(3),zeros(3,9)]; % X5 = u5
X5 = Segment(4).Informed.Q(1:3,1,:); % u5 = X5 of SCS
% Y5
inB5 = inv(Segment(5).Informed.B);
NY5 = [inB5(1,2)*eye(3),inB5(2,2)*eye(3),-inB5(2,2)*eye(3),inB5(3,2)*eye(3)];
Y5 = Mprod_array3(repmat(NY5,[1,1,n]),Segment(5).Informed.Q); % Y5 of SCS
% Z5
NZ5 = [inB5(1,3)*eye(3),inB5(2,3)*eye(3),-inB5(2,3)*eye(3),inB5(3,3)*eye(3)];
Z5 = Mprod_array3(repmat(NZ5,[1,1,n]),Segment(5).Informed.Q); % Z5 of SCS
%
Joint(5).Informed.Kk(1,37:48,:) = [zeros(1,3,n),permute(-X5,[2,1,3]),zeros(1,6,n)] + ...
    Mprod_array3(permute(rV16 - Segment(5).Informed.Q(4:6,1,:),[2,1,3]),repmat(NX5,[1,1,n]));
Joint(5).Informed.Kk(2,37:48,:) = [zeros(1,3,n),permute(-Y5,[2,1,3]),zeros(1,6,n)] + ...
    Mprod_array3(permute(rV16 - Segment(5).Informed.Q(4:6,1,:),[2,1,3]),repmat(NY5,[1,1,n]));
Joint(5).Informed.Kk(3,37:48,:) = [zeros(1,3,n),permute(-Z5,[2,1,3]),zeros(1,6,n)] + ...
    Mprod_array3(permute(rV16 - Segment(5).Informed.Q(4:6,1,:),[2,1,3]),repmat(NZ5,[1,1,n]));

