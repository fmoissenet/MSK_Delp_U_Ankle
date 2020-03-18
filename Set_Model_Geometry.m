% MAIN PROGRAM
% Set_Model_Geometry.m
%__________________________________________________________________________
%
% PURPOSE
% Set generic model geometry
%
% SYNOPSIS
% (Segment,Model) = Scale_Model(Segment)
%
% INPUT
% Segment (cf. data structure in user guide)
%
% OUTPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% DESCRIPTION
% Set generic model geometry
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
%
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Florent Moissenet
% March 2020
%
%__________________________________________________________________________

function [Segment,Model] = Set_Model_Geometry(Segment)

% -------------------------------------------------------------------------
% Generic model geometry (data from Delp 1990 model - LowerExtremityModel)
% -------------------------------------------------------------------------
Model.geometry.height = 1.80;                                 % The Delp model represents a subject that is about 1.8m
Model.geometry.Wm(2) = 0.1070;                                % Foot width between metatarses (m1/m5)
Model.geometry.Lm(2) = 0.1251;                                % Foot length defined from mid-malleolus to mid-metatarses (m1/m5)
Model.geometry.Wm(3) = 0.0761;                                % Width in meter of shank segment
Model.geometry.Lm(3) = 0.3914;                                % Length in meter of shank segment (shank length is 0.3907 but mid-epicondyle to mid-malleolus is -0.3914 on Y-axis)
% Model.geometry.Wm(4) not defined
% Model.geometry.Lm(4) not defined
Model.geometry.Wm(5) = 0.0786;                                % Width in meter of thigh segment
Model.geometry.Lm(5) = 0.4020;                                % Length in meter of thigh segment (thigh length (same as tx + ty at 0° of knee flexion) is 0.3960 but mid-epicondyles is -0.4020 on Y-axis)
% Model.geometry.Wm(6) not defined
% Model.geometry.Lm(6) not defined
Model.geometry.T_calca_talus = [-0.04877; -0.04195; 0.00792]; % Translation (tx, ty, tz) from calcaneus model SCS centre to talus model SCS centre
Model.geometry.T_talus_tibia = [0; -0.426; 0];                % Translation (tx, ty, tz) from talus model SCS centre to tibia model SCS centre (-0.43 in delpjnt.jnt but -0.426 in dynamic_light.jnt)
Model.geometry.T_tibia_mmall = [0.0006; 0.3999; -0.0092];     % Translation (tx, ty, tz) from tibia model SCS centre to mid malleolii (rD3 or rP2)
Model.geometry.T_tibia_mepic = [0; 0.0060; 0.0005];           % Translation (tx, ty, tz) from tibia model SCS centre to mid epicondyles (rD5 or rP3)
Model.geometry.T_pelvi_femur = [0.0707; 0.0661; -0.0835];     % Translation (tx, ty, tz) from pelvis model SCS centre to femur model SCS centre (rD6 or rP5)
Model.geometry.T_patel_patel = [0; -(0.0264+0.025/2); 0];     % Translation (tx, ty, tz) from patella rD4 to half patella heigth on Y-axis
Model.geometry.L_tibialTuber = [0.039; -0.082; 0.000];        % Landmark for tibial tuberosity
