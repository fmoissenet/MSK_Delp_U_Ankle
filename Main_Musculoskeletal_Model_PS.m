% MAIN PROGRAM
% Main_Musculoskeletal_Model_PS.m
%__________________________________________________________________________
%
% PURPOSE
% Building of lower limb model and computation of musculo-tendon, contact,
% ligament and bone forces
%
% SYNOPSIS
% N/A (i.e., main program)
%
% DESCRIPTION
%
% REFERENCE
% F Moissenet, L Cheze, R Dumas. A 3D lower limb musculoskeletal model for
% simultaneous estimation of musculo-tendon, joint contact, ligament and
% bone forces during gait. Journal of Biomechanics 2014;47(1):50-8.
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Modify_Segment_PS.m
% Multibody_Kinematics_Optimisation_PS.m
% Compute_J.m
% Compute_G.m
% Compute_P.m
% Compute_R.m
% Compute_L_Hill.m
% Static_Optimisation_Lagrange_Multipliers_Hill.m
% Main_Muscle_Lines_Visualisation_PS.m
% 
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphael Dumas, Florent Moissenet, Edouard Jouan
% September 2012
%
% Modified by Raphael Dumas
% April 2013
% Computation of the non-optimized Lagrange multipiers
% Updated figure
%
% Modified by Raphael Dumas
% February 2015
% Jacobian of the kinematic constraints for spherical joints in SCS and
% full selection of lambda1
%
% Modified by Raphael Dumas
% January 2018
% Renamed: Multibody Kinematics Optimisation
%
% Modified by Raphael Dumas
% June 2019
% Change w axis into Z axis for thigh segment
% Test and save contact forces
%
% Modified by Raphael Dumas
% October 2019
% Functions _GC5
%
% Modified by Raphael Dumas
% January 2020
% Renamed _PS
%__________________________________________________________________________


% -------------------------------------------------------------------------
% LOAD RECORDED DATA
% -------------------------------------------------------------------------

% Clear/clean workspace, command window and figures
clc
% close all
clear all

% Select data file
[a,b] = uigetfile('*.mat','Selet data file');
load([b,a]);

% Number of frames
n = size(Segment(2).rM,3);

% Informed segment length
Segment(3).L = 0.4316; % From PS CT-scan
Segment(5).L = 0.4256; % From PS CT-scan
 
% Informed segment width
% Segment(3).W = 0.0695; % From PS CT-scan
Segment(5).W = 0.0871; % From PS CT-scan

% Modify point of expression of Joint(1).M
CoP  = Segment(1).Q(4:6,:,:);
CoPs = Segment(2).Q(4:6,:,:); % CoP_star: new point of expression of Joint(1).M
Segment(1).Q(4:6,:,:) = CoPs;
Joint(1).M = Joint(1).M + cross((CoP-CoPs),Joint(1).F);

% -------------------------------------------------------------------------
% COMPUTE POSITIONS, ACCELERATIONS AND JACOBIAN MATRIX
% -------------------------------------------------------------------------

% Insert patella as segment 4
Segment = Modify_Segment_PS(Segment);

% Change w axis into Z axis for thigh segment
Segment(5).T = Q2Tuv_array3(Segment(5).Q);
Segment(5).Q(10:12,:,:) = Segment(5).T(1:3,3,:);

% Re-definition of hip joint centre based on Delp
% LASIS = [0.0; 0.0; -0.128];
% RASIS = [0.0; 0.0; 0.128];
% RHJC = [-0.0707; -0.0661; 0.0835];
s = mean(sqrt(sum((Segment(6).rM(:,2,:) - Segment(6).rM(:,1,:)).^2)),3)...
    /(0.128*2); % Scaled by pelvis width
Segment(6).T = Q2Tuv_array3(Segment(6).Q);
rV16 = mean(Segment(6).rM(:,1:2,:),2) + ... % Mid-ASIS
    Mprod_array3(Segment(6).T(1:3,1:3,:), ... 
    repmat([-0.0707; -0.0661; 0.0835]*s,[1,1,n]));
Segment(5).Q(4:6,:,:) = rV16;
Segment(5).Q(1:3,:,:) = Vnorm_array3(cross(rV16 - Segment(5).Q(7:9,:,:), ...
    Segment(5).Q(10:12,:,:)));
Segment(5).Q(10:12,:,:) = Vnorm_array3(cross(Segment(5).Q(1:3,:,:), ...
    rV16 - Segment(5).Q(7:9,:,:)));
Segment(6).Q(7:9,:,:) = rV16 - Segment(6).T(1:3,3,:)*0.0835*s;    

% Multibody kinematics optimisation
[Segment,Joint] = Multibody_Kinematics_Optimisation_PS(Segment,Joint,f);

% Model velocities
Model.dQdt = [Segment(2).dQdt; ... % Foot (12*1*n)
    Segment(3).dQdt; ... % Shank (12*1*n)
    Segment(4).dQdt; ... % Patella (12*1*n)
    Segment(5).dQdt]; % Thigh (12*1*n)

% Model accelerations
Model.d2Qdt2 = [Segment(2).d2Qdt2; ... % Foot (12*1*n)
    Segment(3).d2Qdt2; ... % Shank (12*1*n)
    Segment(4).d2Qdt2; ... % Patella (12*1*n)
    Segment(5).d2Qdt2]; % Thigh (12*1*n)

% Model Jacobian
Model.K = [Joint(2).Kk(:,1:48,:); ... % Ankle (4*48*n)
    Joint(3).Kk(:,1:48,:); ... % Tibio-femoral (5*48*n)
    Joint(4).Kk(:,1:48,:); ... % Patello-femoral (6*48*n)
    Joint(5).Kk(:,1:48,:); ... % Hip (3*48*n)
    Segment(2).Kr(:,1:48,:); ... % Foot (6*48*n)
    Segment(3).Kr(:,1:48,:); ... % Shank (6*48*n)
    Segment(4).Kr(:,1:48,:); ... % Patella (6*48*n)
    Segment(5).Kr(:,1:48,:)]; % Thigh (6*48*n)

% -------------------------------------------------------------------------
% PREPARE SEGMENT KINETICS
% -------------------------------------------------------------------------
Segment = Compute_J(Segment); % Pseudo inertia matrix
[Segment,Model] = Compute_G(Segment,Model); % Mass matrix

% -------------------------------------------------------------------------
% COMPUTE EXTERNAL FORCES
% -------------------------------------------------------------------------
Model = Compute_P(Segment,Model); % Weight
Model = Compute_R(Segment,Joint,Model); % Ground reaction forces

% -------------------------------------------------------------------------
% COMPUTE MUSCULAR LEVER ARMS
% -------------------------------------------------------------------------
[Segment,Model] = Compute_L_Hill(Segment,Model,height/100,2);

% -------------------------------------------------------------------------
% ESTIMATE MUSCULO-TENDON, CONTACT, LIGAMENT AND BONE FORCES USING A 
% CONSTRAINED STATIC OPTIMISATION METHOD
% -------------------------------------------------------------------------

% Optimisation
Model = Static_Optimisation_Lagrange_Multipliers_Hill(Segment,Joint,Model,f);

% Musculo-tendon forces (tensile forces > 0)
Model.Fm = Model.X(1:43,1,1:n); % About lines of action

% Contact forces (reaction forces > 0 i.e., acting on the proximal segment)
Model.Fc = [Model.X(43+1:43+3,1,1:n); ... % 3D ankle contact in foot SCS
    Model.X(43+4:43+8,1,1:n); ... % 3D tibio-femoral medial and lateral (not on Z) contact in shank SCS
    Model.X(43+9:43+11,1,1:n); ... % 3D patello-femoral contact in patella SCS
    Model.X(43+13:43+15,1,1:n)]; % 3D hip contact in thigh SCS

% Ligament forces (tensile forces > 0)
% About lines of action
Model.Fl = Model.X(43+12,1,1:n); % PT

% Bone forces (compression forces > 0)
% About segment Y axis
Model.Fb = Model.X(43+16:43+19,1,1:n); % Foot, tibia, patella, femur axial

% Contributions
% PS_ngait_og_ss1: 2:84, PS_ngait_og_ss8: 1:84, 
% PS_ngait_og_ss9: 4:83, PS_ngait_og_ss11: 5:84
tTime        = 5:84; % Range of time (stance only) for proper CoP data 
w            = 1e-4; % Weight parameter used in the massed pseudo-inverse
Contribution = computeContributions_Delp_U_Ankle(Segment,Joint,Model,w,tTime,CoP,CoPs);
save([regexprep(a,'.mat',''),'_CoP.mat'],'Contribution');

% -------------------------------------------------------------------------
% FIGURES
% -------------------------------------------------------------------------

% Figure 1
Main_Joint_Kinematics
Main_Muscle_Lines_Visualisation_PS

% Figure 2
figure
plot(squeeze(1/(mass*9.81)*Model.Fm(:,:,:))');
legend({'Gluteus maximus I', 'Gluteus maximus II', 'Gluteus maximus III', ...
    'Gluteus medius I', 'Gluteus medius II', 'Gluteus medius III', ...
    'Gluteus minimus I', 'Gluteus minimus II', 'Gluteus minimus III', ...
    'Adductor longus', 'Adductor brevis', ...
    'Adductor magnus I', 'Adductor magnus II', 'Adductor magnus III', ...
    'Pectineus', 'Illiacus', 'Psoas', ...
    'Quadratus femoris', 'Gemelli', ...
    'Piriformis','Tensor fasciae latae', ...
    'Gracilis', 'Sartorius', 'Semimembranosus', ...
    'Semitendinus', 'Biceps femoris long head', ...
    'Biceps femoris short head', 'Rectus femoris', ...
    'Vastus medialis', 'Vastus intermedialis', 'Vastus lateralis', ...
    'Gastrocnemius medialis', 'Gastrocnemius lateralis', ...
    'Soleus', 'Tibialis posterior', 'Tibialis anterior', ...
    'Peroneus brevis', 'Peroneus longus', 'Peroneus tertius', ...
    'Extensor digitorum longus', 'Extensor hallucis longus', ...
    'Flexor digitorum longus', 'Flexor hallucis longus'})

% Figure 2
figure
plot(squeeze(1/(mass*9.81)*[Model.Fc(:,:,:); ...
    Model.Fl(:,:,:); ...
    Model.Fb(:,:,:)])');
legend({'Ankle contact (about X axis of Foot)', ...
    'Ankle contact (about Y axis of Foot)', ...
    'Ankle contact (about Z axis of Foot)',...
    'Tibio-femoral medial contact (about X axis of Shank)', ...
    'Tibio-femoral medial contact (about Y axis of Shank)', ...
    'Tibio-femoral medial contact (about Z axis of Shank)', ...
    'Tibio-femoral lateral contact (about X axis of Shank)', ...
    'Tibio-femoral lateral contact (about Y axis of Shank)', ...
    'Patello-femoral contact (about X axis of Patella)', ...
    'Patello-femoral contact (about Y axis of Patella)', ...
    'Patello-femoral contact (about Z axis of Patella)', ...
    'Hip contact (about X axis of Thigh)', ...
    'Hip contact (about Y axis of Thigh)', ...
    'Hip contact (about Z axis of Thigh)', ...
    'PT', ...
    'Foot axial', 'Tibia axial', 'Patella axial', 'Femur axial'})

