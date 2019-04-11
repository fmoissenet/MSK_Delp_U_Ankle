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
% Modify_Segment.m
% Multibody_Kinematics_Optimisation.m
% Compute_J.m
% Compute_G.m
% Compute_P.m
% Compute_R.m
% Compute_L.m
% Static_Optimisation_Lagrange_Multipliers.m
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
% Jacobian of the kinematic constraints for spherical joints in SCS
% and full selection of lambda1
%
% Modified by Raphael Dumas
% August 2017
% Universal joint at the ankle
%
% Modified by Raphael Dumas
% January 2018
% Renamed: Multibody Kinematics Optimisation
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

addpath('C:\Users\florent.moissenet\Documents\Professionnel\routines\github\MSK_Delp_U_Ankle');
cd('C:\Users\florent.moissenet\Documents\Professionnel\routines\github\MSK_Delp_U_Ankle');

% -------------------------------------------------------------------------
% LOAD RECORDED DATA
% -------------------------------------------------------------------------

% Clear/clean workspace, command window and figures
clc
close all
clear all

% Select data file
load('C:\Users\florent.moissenet\Documents\Professionnel\routines\github\MSK_Delp_U_Ankle\example\gait2.mat');
mass = weight;
% Segment = Multibody_Optimisation_Lower_Limb(Segment,'SSS');
% Segment = Multibody_Optimisation_Lower_Limb(Segment,'UnSS');
Segment = Multibody_Optimisation_Lower_Limb(Segment,'UnPS');

% Number of frames and frequency
n = size(Segment(2).rM,3);
f = freq;

% -------------------------------------------------------------------------
% COMPUTE POSITIONS, ACCELERATIONS AND JACOBIAN MATRIX
% -------------------------------------------------------------------------

% Insert patella as segment 4
Segment = Modify_Segment(Segment);
% Optimisation
[Segment,Joint] = Multibody_Kinematics_Optimisation(Segment,Joint,f);

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
[Segment,Model] = Compute_L(Segment,Model);

% -------------------------------------------------------------------------
% ESTIMATE MUSCULO-TENDON, CONTACT, LIGAMENT AND BONE FORCES USING A 
% CONSTRAINED STATIC OPTIMISATION METHOD
% -------------------------------------------------------------------------

% Optimisation
Model = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model);

% Musculo-tendon forces (tensile forces > 0)
Model.Fm = Model.X(1:43,1,1:n); % About lines of action

% Contact forces (reaction forces > 0 i.e., acting on the proximal segment)
Model.Fc = [Model.X(43+1:43+3,1,1:n); ... % 3D ankle contact in foot SCS
    Model.X(43+4,1,1:n);Model.X(43+5,1,1:n); ... % Tibio-femoral medial and lateral contact in shank SCS
    Model.X(43+9:43+11,1,1:n); ... % 3D patello-femoral contact in patella SCS
    Model.X(43+13:43+15,1,1:n)]; % 3D hip contact in thigh SCS

% Ligament forces (tensile forces > 0)
% About lines of action
Model.Fl = [Model.X(43+6:43+8,1,1:n); ... % ACL, PCL, MCL
    Model.X(43+12,1,1:n)]; % PT

% Bone forces (compression forces > 0)
% About segment Y axis
Model.Fb = Model.X(43+16:43+19,1,1:n); % Foot, tibia, patella, femur axial

% Contributions
Contribution = computeContributions_Delp_U_Ankle(Segment,Joint,Model,weight);
%%
% -------------------------------------------------------------------------
% FIGURES
% -------------------------------------------------------------------------

% Figure 1
Main_Joint_Kinematics;

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

% Figure 3
% Figure 3
figure
plot(squeeze(1/(mass*9.81)*[Model.Fc(:,:,:); ...
    Model.Fl(:,:,:); ...
    Model.Fb(:,:,:)])');
legend({'Ankle contact (about X axis of Foot)', 'Ankle contact (about Y axis of Foot)', 'Ankle contact (about Z axis of Foot)',...
    'Tibio-femoral medial contact', 'Tibio-femoral lateral contact', ...
    'Patello-femoral contact (about X axis of Patella)', 'Patello-femoral contact (about Y axis of Patella)', 'Patello-femoral contact (about Z axis of Patella)', ...
    'Hip contact (about X axis of Thigh)', 'Hip contact (about Y axis of Thigh)', 'Hip contact (about Z axis of Thigh)', ...
    'ACL', 'PCL', 'MCL', 'PT', ...
    'Foot axial', 'Tibia axial', 'Patella axial', 'Femur axial'})