% MAIN PROGRAM
% Main_Musculoskeletal_Model.m
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
% Insert_Patella.m
% Multibody_Kinematics_Optimisation.m
% Compute_J.m
% Compute_G.m
% Compute_P.m
% Compute_R.m
% Compute_L_Hill.m
% Static_Optimisation_Lagrange_Multipliers_Hill.m
% Main_Muscle_Lines_Visualisation.m
% 
% MATLAB VERSION
% Matlab R2020a
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
% Renamed to multibody kinematics optimisation
%
% Modified by Raphael Dumas
% June 2019
% Change w axis into Z axis for thigh segment
%
% Modified by Raphael Dumas
% March 2020
% Generic vs. Informed structures (Segment, Joint, Model)
% n,f,fc as Model.Informed fields
% Modify_Segment.m renamed Insert_Patella.m
%__________________________________________________________________________


%% -------------------------------------------------------------------------
% LOAD RECORDED DATA
% -------------------------------------------------------------------------

% Clear/clean workspace, command window and figures
clc
close all
clear all

% Select data file
[a,b] = uigetfile('*.mat','Selet data file');
load([b,a]);

% Antropometry
try
    Model.Informed.Mass = mass;
catch
    Model.Informed.Mass = weight; % Wrong variable name
end
Model.Informed.Height = height/100; % in m

% Number of frames
Model.Informed.n = n;
% Acquisition frequency
Model.Informed.f = f;
% Cut off frequency for filtering
Model.Informed.fc = 5;

% Informed segment parameters
for i = 1:3 % From i = 1 (Forceplate) to i = 3 (Shank)
    Segment(i).Informed = Segment(i); % Add all fields at a second level
end
% Prepare insertion for patella
Segment(5).Informed = Segment(4); % Add all fields at a second level
Segment(6).Informed = Segment(5); % Add all fields at a second level
%
Segment(2).Informed = rmfield(Segment(2).Informed,'Informed'); 
Segment(3).Informed = rmfield(Segment(3).Informed,'Informed');
Segment(5).Informed = rmfield(Segment(5).Informed,'Informed');
Segment(6).Informed = rmfield(Segment(6).Informed,'Informed');
%
% i = 1 (Ground)
Joint(1).Informed = Joint(1); % Add all fields at a second level

% Original fields are not kept
Segment = rmfield(Segment,'Q');
Segment = rmfield(Segment,'rM');
Segment = rmfield(Segment,'m');
Segment = rmfield(Segment,'rCs');
Segment = rmfield(Segment,'Is');
Joint = rmfield(Joint,'F');
Joint = rmfield(Joint,'M');


%% -------------------------------------------------------------------------
% SET GEOMETRY
% -------------------------------------------------------------------------

% Defined according to Delp's model
Main_Set_Geometry_Delp

% Optional
% Change lumbar and hip joint centres to better match Delp's one
% Based on real markers rM16 = LASIS and rM26 = ASIS
Segment(6).Informed.W = mean(sqrt(sum((Segment(6).Informed.rM(:,1,:) - Segment(6).Informed.rM(:,2,:)).^2)),3);
Segment(6).Informed.Scale = Segment(6).Informed.W/...
    Segment(6).Generic.W; % Default scale
Segment(6).Informed.T = Q2Tuv_array3(Segment(6).Informed.Q);
Segment(6).Informed.Q(4:6,:,:) = ...
    (Segment(6).Informed.rM(:,1,:) + Segment(6).Informed.rM(:,2,:))/2 - ... % Mid-ASIS
    0.01*Segment(6).Informed.T(1:3,1,:) + ... Flesh margin on X-axis (1 cm)
    Mprod_array3(Segment(6).Informed.T(1:3,1:3,:), ... 
    repmat([-t_lumbar2pelvis]*Segment(6).Informed.Scale,[1,1,n])); % Redefined P6
Segment(6).Informed.Q(7:9,:,:) = Segment(6).Informed.Q(4:6,:,:) + ... 
    Mprod_array3(Segment(6).Informed.T(1:3,1:3,:), ... 
    repmat([Segment(6).Generic.rVs(1:2,1);0]*Segment(6).Informed.Scale,[1,1,n])); % Redefined D6
Segment(5).Informed.Q(4:6,:,:) = Segment(6).Informed.Q(7:9,:,:) + ...
    Mprod_array3(Segment(6).Informed.T(1:3,1:3,:), ... 
    repmat([0;0;Segment(6).Generic.rVs(3,1)]*Segment(6).Informed.Scale,[1,1,n])); % Redefined P5
Segment(5).Informed.Q(1:3,:,:) = Vnorm_array3(cross(Segment(5).Informed.Q(4:6,:,:) ...
    - Segment(5).Informed.Q(7:9,:,:), ...
    Segment(5).Informed.Q(10:12,:,:)));  % Rededined u5 axis
Segment(5).Informed.Q(10:12,:,:) = Vnorm_array3(cross(Segment(5).Informed.Q(1:3,:,:), ...
    Segment(5).Informed.Q(4:6,:,:) - Segment(5).Informed.Q(7:9,:,:))); % Rededined w5 axis

% Define segment scales
[Segment,Model] = Define_Scale(Segment,Model);


%% -------------------------------------------------------------------------
% COMPUTE POSITIONS, ACCELERATIONS AND JACOBIAN MATRIX
% -------------------------------------------------------------------------

% Define patella
Segment = Insert_Patella(Segment,Model);

% Multibody kinematics optimisation
[Segment,Joint] = Multibody_Kinematics_Optimisation(Segment,Joint,Model);
Main_Joint_Kinematics

% Model velocities
Model.Informed.dQdt = [Segment(2).Informed.dQdt; ... % Foot (12*1*n)
    Segment(3).Informed.dQdt; ... % Shank (12*1*n)
    Segment(4).Informed.dQdt; ... % Patella (12*1*n)
    Segment(5).Informed.dQdt]; % Thigh (12*1*n)

% Model accelerations
Model.Informed.d2Qdt2 = [Segment(2).Informed.d2Qdt2; ... % Foot (12*1*n)
    Segment(3).Informed.d2Qdt2; ... % Shank (12*1*n)
    Segment(4).Informed.d2Qdt2; ... % Patella (12*1*n)
    Segment(5).Informed.d2Qdt2]; % Thigh (12*1*n)

% Model Jacobian
Model.Informed.K = [Joint(2).Informed.Kk(:,1:48,:); ... % Ankle (4*48*n)
    Joint(3).Informed.Kk(:,1:48,:); ... % Tibio-femoral (5*48*n)
    Joint(4).Informed.Kk(:,1:48,:); ... % Patello-femoral (6*48*n)
    Joint(5).Informed.Kk(:,1:48,:); ... % Hip (3*48*n)
    Segment(2).Informed.Kr(:,1:48,:); ... % Foot (6*48*n)
    Segment(3).Informed.Kr(:,1:48,:); ... % Shank (6*48*n)
    Segment(4).Informed.Kr(:,1:48,:); ... % Patella (6*48*n)
    Segment(5).Informed.Kr(:,1:48,:)]; % Thigh (6*48*n)


%% -------------------------------------------------------------------------
% PREPARE SEGMENT KINETICS
% -------------------------------------------------------------------------
Segment = Compute_J(Segment,Model); % Pseudo inertia matrix
[Segment,Model] = Compute_G(Segment,Model); % Mass matrix


%% -------------------------------------------------------------------------
% COMPUTE EXTERNAL FORCES
% -------------------------------------------------------------------------
Model = Compute_P(Segment,Model); % Weight
Model = Compute_R(Segment,Joint,Model); % Ground reaction forces


%% -------------------------------------------------------------------------
% COMPUTE MUSCULAR LEVER ARMS
% -------------------------------------------------------------------------
[Segment,Model] = Compute_L_Hill(Segment,Model);
Main_Muscle_Lines_Visualisation


%% -------------------------------------------------------------------------
% ESTIMATE MUSCULO-TENDON, CONTACT, LIGAMENT AND BONE FORCES USING A 
% CONSTRAINED STATIC OPTIMISATION METHOD
% -------------------------------------------------------------------------

% Optimisation
Model = Static_Optimisation_Lagrange_Multipliers_Hill(Segment,Joint,Model);

% Musculo-tendon forces (tensile forces > 0)
Model.Informed.Fm = Model.Informed.X(1:43,1,1:n); % About lines of action
% Muscle activation (0 < a < 1)
Model.Informed.a = Mprod_array3(Model.Informed.A,Model.Informed.Fm) + ...
    Model.Informed.B;

% Contact forces (reaction forces > 0 i.e., acting on the proximal segment)
Model.Informed.Fc = [Model.Informed.X(43+1:43+3,1,1:n); ... % 3D ankle contact in foot SCS
    Model.Informed.X(43+4:43+8,1,1:n); ... % 3D tibio-femoral medial and lateral (not on Z) contact in shank SCS
    Model.Informed.X(43+9:43+11,1,1:n); ... % 3D patello-femoral contact in patella SCS
    Model.Informed.X(43+13:43+15,1,1:n)]; % 3D hip contact in thigh SCS

% Ligament forces (tensile forces > 0)
% About lines of action
Model.Informed.Fl = Model.Informed.X(43+12,1,1:n); % PT

% Bone forces (compression forces > 0)
% About segment Y axis
Model.Informed.Fb = Model.Informed.X(43+16:43+19,1,1:n); % Foot, tibia, patella, femur axial


%% -------------------------------------------------------------------------
% FIGURES
% -------------------------------------------------------------------------

% Musculo-tendon forces
figure
plot(squeeze(1/(Model.Informed.Mass*9.81)*Model.Informed.Fm(:,:,:))');
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

figure
plot(squeeze(Model.Informed.a(:,:,:))');
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

% Contact and ligament forces 
figure
plot(squeeze(1/(Model.Informed.Mass*9.81)*...
    [Model.Informed.Fc(:,:,:); ...
    Model.Informed.Fl(:,:,:); ...
    Model.Informed.Fb(:,:,:)])');
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

