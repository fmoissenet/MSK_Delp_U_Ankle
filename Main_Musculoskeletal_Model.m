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


% -------------------------------------------------------------------------
% LOAD RECORDED DATA
% -------------------------------------------------------------------------

% Clear/clean workspace, command window and figures
clc;
close all;
clearvars;
cd('C:\Users\florent.moissenet\Documents\Professionnel\routines\github\MSK_Delp_U_Ankle');
addpath('C:\Users\florent.moissenet\Documents\Professionnel\routines\github\MSK_Delp_U_Ankle');
addpath('C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Multi-objective optimisation\data\');

grandChallenge = 1;
filename = {'jw_ngait_2' 'jw_ngait_3' 'jw_ngait_4' 'jw_ngait_5' 'jw_ngait_6'};
% grandChallenge = 2;
% filename = {'dm_ngait4' 'dm_ngait10' 'dm_ngait11' 'dm_ngait12' 'dm_ngait13'};
% grandChallenge = 3;
% filename = {'SC_ngait_og5' 'SC_ngait_og6' 'SC_ngait_og7' 'SC_ngait_og8' 'SC_ngait_og9'};
% grandChallenge = 5;
% filename = {'PS_ngait_og_ss1' 'PS_ngait_og_ss7' 'PS_ngait_og_ss8' 'PS_ngait_og_ss9' 'PS_ngait_og_ss11'};

for i = 1:length(filename)
clearvars -except grandChallenge filename i weights;
load(['C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Multi-objective optimisation\data\grand_challenge_',num2str(grandChallenge),'\',filename{i},'_pPS.mat']);

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

% % Optimisation - Weight sum method 1
% weights = zeros(19,1);
% Model0 = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model,Force,weights,'wsm');
% Model0.Fm = Model0.X(1:43,1,1:n); % About lines of action
% Model0.Fc = [Model0.X(43+1:43+3,1,1:n); ... % 3D ankle contact in foot SCS
%     Model0.X(43+4,1,1:n);Model0.X(43+5,1,1:n); ... % Tibio-femoral medial and lateral contact in shank SCS
%     Model0.X(43+9:43+11,1,1:n); ... % 3D patello-femoral contact in patella SCS
%     Model0.X(43+13:43+15,1,1:n)]; % 3D hip contact in thigh SCS
% Model0.Fl = [Model0.X(43+6:43+8,1,1:n); ... % ACL, PCL, MCL
%     Model0.X(43+12,1,1:n)]; % PT
% Model0.Fb = Model0.X(43+16:43+19,1,1:n); % Foot, tibia, patella, femur axial
% Model0 = Dickerson_new(Model0,Emg,'Delp',n);
% 
% % Optimisation - Weight sum method 2
% weights = [1e0;1e0;1e0;2e0;4e0;1e-6;1e-6;1e-6;1e0;1e0;1e0;1e-6;1e-6;1e-6;1e-6;1e-6;1e-6;1e-6;1e-6];
% Model1 = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model,Force,weights,'wsm');
% Model1.Fm = Model1.X(1:43,1,1:n); % About lines of action
% Model1.Fc = [Model1.X(43+1:43+3,1,1:n); ... % 3D ankle contact in foot SCS
%     Model1.X(43+4,1,1:n);Model1.X(43+5,1,1:n); ... % Tibio-femoral medial and lateral contact in shank SCS
%     Model1.X(43+9:43+11,1,1:n); ... % 3D patello-femoral contact in patella SCS
%     Model1.X(43+13:43+15,1,1:n)]; % 3D hip contact in thigh SCS
% Model1.Fl = [Model1.X(43+6:43+8,1,1:n); ... % ACL, PCL, MCL
%     Model1.X(43+12,1,1:n)]; % PT
% Model1.Fb = Model1.X(43+16:43+19,1,1:n); % Foot, tibia, patella, femur axial
% Model1 = Dickerson_new(Model1,Emg,'Delp',n);

% Optimisation - Min max method
weights = zeros(19,1);
Model2 = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model,Force,weights,'mmm');
Model2.Fm = Model2.X(1:43,1,1:n); % About lines of action
Model2.Fc = [Model2.X(43+1:43+3,1,1:n); ... % 3D ankle contact in foot SCS
    Model2.X(43+4,1,1:n);Model2.X(43+5,1,1:n); ... % Tibio-femoral medial and lateral contact in shank SCS
    Model2.X(43+9:43+11,1,1:n); ... % 3D patello-femoral contact in patella SCS
    Model2.X(43+13:43+15,1,1:n)]; % 3D hip contact in thigh SCS
Model2.Fl = [Model2.X(43+6:43+8,1,1:n); ... % ACL, PCL, MCL
    Model2.X(43+12,1,1:n)]; % PT
Model2.Fb = Model2.X(43+16:43+19,1,1:n); % Foot, tibia, patella, femur axial
Model2 = Dickerson_new(Model2,Emg,'Delp',n);

save(['C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Multi-objective optimisation\data\grand_challenge_',num2str(grandChallenge),'\',filename{i},'_results.mat']);

end

outcomes = processing(grandChallenge,filename);
save(['C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Multi-objective optimisation\data\grand_challenge_',num2str(grandChallenge),'\goodnessOfFit_results.mat']);