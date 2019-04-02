% FUNCTION
% Multibody_Optimisation_Lower_Limb.m
%__________________________________________________________________________
%
% PURPOSE
% Multibody optimisation : minimisation of the sum of the squared distances
% between measured and model-determined marker positions subject to
% kinematic and rigid body constraints
%
% SYNOPSIS
% Segment = Multibody_Optimisation_Lower_Limb(Segment,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model (type of kinematic constraints)
%
% OUTPUT
% Segment (cf. data structure in user guide)
%
% DESCRIPTION
% Computation of Q by minimisation under constraints (by Gauss-Newton)
%
% Joint models are:
% 'N': None (6 degrees of freedom)
% 'S': Spherical (3 degrees of freedom)
% 'Un': Universal with neutral alignment (2 degrees of freedom), not for hip
% 'Um': Universal with mean alignment (2 degrees of freedom), not for hip
% 'Hn': Hinge with neutral alignment (1 degree of freedom), not for ankle and hip
% 'Hm': Hinge with mean alignment (1 degree of freedom), not for ankle and hip
% 'P': Parallel mechanism (i.e., Deltazero in Gasparutto et al. 2015), not for hip
% 'Lm': Minimized ligament length variations (i.e., Deltamin in Gasparutto et al. 2015), not for ankle and hip
% 'Ln': Targeted (nominal) ligament length variations (i.e., Deltatheta in Gasparutto et al. 2015), not for ankle and hip
%
% REFERENCES
% S Duprey, L Cheze, R Dumas. Influence of joint constraints on lower
% limb kinematics estimation from skin markers using global optimization.
% Journal of Biomechanics 2010;43(14):2858-2862.
% X Gasparutto, N Sancisi, E Jacquelin, V Parenti-Castelli, R Dumas.
% Validation of a multi-body optimization with knee kinematic models 
% including ligament constraints. Journal of Biomechanics 2015; 48(6):1141
% -1146
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX)
% Vnop_array3.m
% Mprod_array3.m
% Minv_array3.m
%
% MATLAB VERSION
% Matlab R2013a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas
% August 2013
%
% Modified by Raphaël Dumas
% April 2015
% Updated comments
%__________________________________________________________________________

function Segment = Multibody_Optimisation_Lower_Limb(Segment,Model)

% Number of frames
n = size(Segment(2).rM,3);
% Initialisation
Joint(1).Kk = [];


%% ------------------------------------------------------------------------
% Model parameters
% -------------------------------------------------------------------------

% Initialisation
Segment(1).L = NaN; % No value for segment 1 (Forceplate)
Segment(1).alpha = NaN; % No value for segment 1 (Forceplate)
Segment(1).beta = NaN; % No value for segment 1 (Forceplate)
Segment(1).gamma = NaN; % No value for segment 1 (Forceplate)
Joint(1).d(1,1:6) = NaN(1,6); % No value for joint 1 (CoP)

% Mean segment geometry and markers coordinates
for i = 2:5 % From i = 2 (foot) to i = 5 (pelvis)
    
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


%% ------------------------------------------------------------------------
% Joint parameters and initial guess
% -------------------------------------------------------------------------

switch Model % Ankle
    
    case {'NNN','NNS','NHnN','NHnS','NHmN','NHmS','NUnN','NUnS','NUmN', ...
            'NUmS','NSN','NSS','NPN','NPS',...
            'NLmN','NLnN','NLmS','NLnS'} % No ankle joint
        
        % Initial guess for Lagrange multipliers
        lambdakA = []; % To be concatenated
        
    case {'UnNN','UnNS','UnHnN','UnHnS','UnHmN','UnHmS','UnUnN', ...
            'UnUnS','UnUmN','UnUmS','UnSN','UnSS','UnPN','UnPS',...
            'UnLmN','UnLnN','UnLmS','UnLnS'} % Universal ankle joint
        
        % Ankle U-joint angle from neutral position
        thetaA = Segment(3).beta; % u2 = u3 at neutral position
        % Initial guess for Lagrange multipliers
        lambdakA = zeros(4,1,n);
        
    case {'UmNN','UmNS','UmHnN','UmHnS','UmHmN','UmHmS','UmUnN', ...
            'UmUnS','UmUmN','UmUmS','UmSN','UmSS','UmPN','UmPS',...
            'UmLmN','UmLnN','UmLmS','UmLnS'} % Universal ankle joint
        
        % Ankle U-joint mean angle
        thetaA = mean(acos(dot(Segment(3).Q(10:12,1,:), ...
            Segment(2).Q(1:3,1,:))),3)*180/pi;
        % Initial guess for Lagrange multipliers
        lambdakA = zeros(4,1,n);
        
    case {'SNN','SNS','SHnN','SHnS','SHmN','SHmS','SUnN','SUnS','SUmN', ...
            'SUmS','SSN','SSS','SPN','SPS',...
            'SLmN','SLnN','SLmS','SLnS'} % Spherical ankle joint
        
        % Initial guess for Lagrange multipliers
        lambdakA = zeros(3,1,n);
        
    case {'PNN','PNS','PHnN','PHnS','PHmN','PHmS','PUnN','PUnS', ...
            'PUmN','PUmS','PSN','PSS','PPN','PPS', ...
            'PLmN','PLmS','PLnN','PLnS'} % Paralell ankle mechanism
        
        % Generic model of the right ankle
        % Di Gregorio et al. Med Bio Eng Comput 2007, 45: 305-313
        % Mean of specimen #1 to #3 (with mirror symetry for left specimens)
        % In talus/calcaneus SCS
        % Origin at midpoint between medial and lateral centers in calcaneus
        % Second proposed parallel mechanism (1 spherical joint & 2 ligaments)
        
        % Virtual markers
        rVs12 = [0.001553; 0.000275; 0.009531]; % Lateral_centre_calcaneus = Lateral_contact_tibia
        rVs22 = [0.001653; -0.017121; -0.015714]; % TiCaL_calcaneus
        rVs32 = [-0.022070; -0.016800; 0.025078]; % CaFiL_calcaneus
        rVs13 = rVs12; % Lateral_contact_tibia
        rVs23 = [-0.000007; 0.016023; -0.026644]; % TiCaL_tibia
        rVs33 = [0.001954; -0.001809; 0.027941]; % CaFiL_fibula
        
        % Ligament lengths: distances (associated with joint 2)
        Joint(2).d(1,1) = norm(rVs22 - rVs23); % TiCaL length
        Joint(2).d(1,2) = norm(rVs32 - rVs33); % CaFiL length
        
        % Coordinates of virtual markers in segment 2
        % Expressed in (u2, rP2-rD2, w2)
        Segment(2).nV(:,1) = inv(Segment(2).B)*rVs12; % V12: virtual marker 1 of segment 2
        Segment(2).nV(:,2) = inv(Segment(2).B)*rVs22; % V22: virtual marker 2 of segment 2
        Segment(2).nV(:,3) = inv(Segment(2).B)*rVs32; % V32: virtual marker 3 of segment 2
        
        % Coordinates of virtual markers in segment 3
        % Expressed in (u3, rP3-rD3, w3)
        Segment(3).nV(:,1) = [0;-1;0] + inv(Segment(3).B)*rVs13; % V13: virtual marker 1 of segment 3
        Segment(3).nV(:,2) = [0;-1;0] + inv(Segment(3).B)*rVs23; % V23: virtual marker 2 of segment 3
        Segment(3).nV(:,3) = [0;-1;0] + inv(Segment(3).B)*rVs33; % V33: virtual marker 3 of segment 3
        
        % Interpolation matrices
        NV12 = [Segment(2).nV(1,1)*eye(3),...
            (1 + Segment(2).nV(2,1))*eye(3), ...
            - Segment(2).nV(2,1)*eye(3), ...
            Segment(2).nV(3,1)*eye(3)];
        NV22 = [Segment(2).nV(1,2)*eye(3),...
            (1 + Segment(2).nV(2,2))*eye(3), ...
            - Segment(2).nV(2,2)*eye(3), ...
            Segment(2).nV(3,2)*eye(3)];
        NV32 = [Segment(2).nV(1,3)*eye(3),...
            (1 + Segment(2).nV(2,3))*eye(3), ...
            - Segment(2).nV(2,3)*eye(3), ...
            Segment(2).nV(3,3)*eye(3)];
        NV13 = [Segment(3).nV(1,1)*eye(3),...
            (1 + Segment(3).nV(2,1))*eye(3), ...
            - Segment(3).nV(2,1)*eye(3), ...
            Segment(3).nV(3,1)*eye(3)];
        NV23 = [Segment(3).nV(1,2)*eye(3),...
            (1 + Segment(3).nV(2,2))*eye(3), ...
            - Segment(3).nV(2,2)*eye(3), ...
            Segment(3).nV(3,2)*eye(3)];
        NV33 = [Segment(3).nV(1,3)*eye(3),...
            (1 + Segment(3).nV(2,3))*eye(3), ...
            - Segment(3).nV(2,3)*eye(3), ...
            Segment(3).nV(3,3)*eye(3)];
        
        % Initial guess for Lagrange multipliers
        lambdakA = zeros(5,1,n);
        
    otherwise
        display ('No appropriate model for the ankle joint')
        
end

% -------------------------------------------------------------------------
switch Model % Knee
    
    case {'NNN','NNS','UnNN','UnNS','UmNN','UmNS','SNN','SNS','PNN','PNS'} % No knee joint
        
        % Initial guess for Lagrange multipliers
        lambdakK = []; % To be concatenated
        
    case {'NHnN','NNnS','UnHnN','UnHnS','UmHnN','UmHnS','SHnN','SHnS', ...
            'PHnN','PHnS'} % Hinge knee joint
        
        % Knee hinge angles from neutral position
        thetaK1 = Segment(4).alpha; % v3 (rP3 - rD3) = v4 (rP4 - rD4) at neutral position
        thetaK2 = Segment(4).beta; % u3 = u4 at neutral position
        % Initial guess for Lagrange multipliers
        lambdakK = zeros(5,1,n);
        
    case {'NHmN','NHmS','UnHmN','UnHmS','UmHmN','UmHmS','SHmN','SHmS', ...
            'PHmN','PHmS'} % Hinge knee joint
        
        % Knee mean hinge angles
        thetaK1 = mean(acos(dot(Segment(4).Q(10:12,1,:), ...
            Segment(3).Q(4:6,1,:) - Segment(3).Q(7:9,1,:))./...
            sqrt(sum((Segment(3).Q(4:6,1,:) - ...
            Segment(3).Q(7:9,1,:)).^2))),3)*180/pi;
        thetaK2 = mean(acos(dot(Segment(4).Q(10:12,1,:), ...
            Segment(3).Q(1:3,1,:))),3)*180/pi;
        % Initial guess for Lagrange multipliers
        lambdakK = zeros(5,1,n);
        
    case {'NUnN','NUnS','UnUnN','UnUnS','UmUnN','UmUnS','SUnN','SUnS', ...
            'PUnN','PUnS'} % Universal knee joint
        
        % Knee U-joint angles from neutral position
        thetaK1 = Segment(4).alpha; % v3 (rP3 - rD3) = v4 (rP4 - rD4) at neutral position
        % Initial guess for Lagrange multipliers
        lambdakK = zeros(4,1,n);
        
    case {'NUmN','NUmS','UnUmN','UnUmS','UmUmN','UmUmS','SUmN','SUmS', ...
            'PUmN','PUmS'} % Universal knee joint
        
        % Knee mean U-joint angles
        thetaK1 = mean(acos(dot(Segment(4).Q(10:12,1,:), ...
            Segment(3).Q(4:6,1,:) - Segment(3).Q(7:9,1,:))./...
            sqrt(sum((Segment(3).Q(4:6,1,:) - ...
            Segment(3).Q(7:9,1,:)).^2))),3)*180/pi;
        % Initial guess for Lagrange multipliers
        lambdakK = zeros(4,1,n);
        
    case {'NSN','NSS','UnSN','UnSS','UmSN','UmSS','SSN','SSS','PSN','PSS'} % Spherical knee joint
        
        % Initial guess for Lagrange multipliers
        lambdakK = zeros(3,1,n);
        
    case {'NPN','NPS','UnPN','UnPS','UmPN','UmPS','SPN','SPS','PPN','PPS', ... % Parallel knee mechanism (zero ligament length variations)
            'NLmN','NLmS','UnLmN','UnLmS','UmLmN','UmLmS','SLmN','SLmS','PLmN','PLmS', ... % Minimized ligament length variations
            'NLnN','NLnS','UnLnN','UnLnS','UmLnN','UmLnS','SLnN','SLnS','PLnN','PLnS'} % Targeted (nominal) ligament length variations
        
        % Generic model of the right knee
        % Gasparutto et al. J Biomech 2015, 48: 1141–1146
        % In tibia SCS = femur SCS at full extension
        Medial_centre_femur = [0.0002; 0.0034; -0.0232];
        Lateral_centre_femur = [-0.0033; 0.0021; 0.0262];
        ACL_femur = [-0.0068; 0.0075; 0.0092];
        PCL_femur = [-0.0027; -0.0011; -0.0022];
        MCL_femur = [0.0028; 0.0058; -0.0476];
        LCL_femur = [0.0033; 0.0023; 0.0362];
        Medial_contact_tibia = [-0.0021; -0.0286; -0.0191];
        Lateral_contact_tibia = [-0.0028; -0.0261; 0.0244];
        n_Medial = [0.0675; 0.9896; -0.1273]; % Medial tibial plateau
        n_Medial = n_Medial/norm(n_Medial);
        n_Lateral = [-0.0881; 0.9942; 0.0617]; % Lateral tibial plateau
        n_Lateral = n_Lateral/norm(n_Lateral);
        ACL_tibia = [0.0128; -0.0261; -0.0001];
        PCL_tibia = [-0.0259; -0.0381; -0.0035];
        MCL_tibia = [0.0021; -0.1171; -0.0058];
        LCL_tibia = [-0.0243; -0.0480; 0.0371];
        % Ligament lengths in neutral position
        d3_0 = norm(ACL_femur - ACL_tibia); % ACL
        d4_0 = norm(PCL_femur - PCL_tibia); % PCL
        d5_0 = norm(MCL_femur - MCL_tibia); % MCL
        d6_0 = norm(LCL_femur - LCL_tibia); % LCL
        
        % Relative ligament lenght
        % Gasparutto et al. J Biomech 2015, 48: 1141–1146
        a = [-2.5e-3, -1.4023e-5, 3.2187e-6, 9.1037e-8, 1.0491e-9, 5.8532e-12, 1.3559e-14; ... ACL
            2.7e-3, -4.2080e-5, -8.5132e-6, -2.4381e-7, -3.0408e-9, -1.8055e-11, -4.1819e-14; ... PCL
            2.3336e-4, 3.2597e-5, 1.3686e-6, 2.4415e-8, 1.7782e-10, 2.2257e-13, -1.7572e-15; ... MCL
            3.5e-3, -1.5564e-4, -1.9254e-5, -5.7922e-7, -7.9173e-9, -5.1875e-11, -1.3239e-13]; ... LCL
            
        % Distances (associated with segment 3)
        Joint(3).d(1,1,:) = repmat(0.0323,[1 1 n]); % Medial radius
        Joint(3).d(1,2,:) = repmat(0.0283,[1 1 n]); % Lateral radius
        
        % i = 3 (leg)
        % Coordinates of virtual markers in segment 3
        % Expressed in  in (u3, rP3-rD3, w3)
        Segment(3).nV(:,4) = inv(Segment(3).B)*Medial_contact_tibia; % V43: virtual marker 4 of segment 3
        Segment(3).nV(:,5) = inv(Segment(3).B)*Lateral_contact_tibia; % V53: virtual marker 5 of segment 3
        Segment(3).nV(:,6) = inv(Segment(3).B)*ACL_tibia; % V63 virtual marker 6 of segment 3
        Segment(3).nV(:,7) = inv(Segment(3).B)*PCL_tibia; % V73 virtual marker 7 of segment 3
        Segment(3).nV(:,8) = inv(Segment(3).B)*MCL_tibia; % V83 virtual marker 8 of segment 3
        Segment(3).nV(:,9) = inv(Segment(3).B)*LCL_tibia; % V93 virtual marker 9 of segment 3
        
        % Normals (associated with segment 3)
        Segment(3).nn(:,1) = inv(Segment(3).B)*n_Medial;
        Segment(3).nn(:,2) = inv(Segment(3).B)*n_Lateral;
        
        % Coordinates of virtual markers in segment 4
        % Expressed in  in (u4, rP4-rD4, w4)
        Segment(4).nV(:,1) = [0;-1;0] + inv(Segment(4).B)*Medial_centre_femur; % V14: virtual marker 1 of segment 4
        Segment(4).nV(:,2) = [0;-1;0] + inv(Segment(4).B)*Lateral_centre_femur; % V24: virtual marker 2 of segment 4
        Segment(4).nV(:,3) = [0;-1;0] + inv(Segment(4).B)*ACL_femur; % V34: virtual marker 3 of segment 4
        Segment(4).nV(:,4) = [0;-1;0] + inv(Segment(4).B)*PCL_femur; % V44: virtual marker 4 of segment 4
        Segment(4).nV(:,5) = [0;-1;0] + inv(Segment(4).B)*MCL_femur; % V54: virtual marker 5 of segment 4
        Segment(4).nV(:,6) = [0;-1;0] + inv(Segment(4).B)*LCL_femur; % V64: virtual marker 6 of segment 4
        
        % Interpolation matrices
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
        NV14 = [Segment(4).nV(1,1)*eye(3),...
            (1 + Segment(4).nV(2,1))*eye(3), ...
            - Segment(4).nV(2,1)*eye(3), ...
            Segment(4).nV(3,1)*eye(3)];
        NV24 = [Segment(4).nV(1,2)*eye(3),...
            (1 + Segment(4).nV(2,2))*eye(3), ...
            - Segment(4).nV(2,2)*eye(3), ...
            Segment(4).nV(3,2)*eye(3)];
        NV34 = [Segment(4).nV(1,3)*eye(3),...
            (1 + Segment(4).nV(2,3))*eye(3), ...
            - Segment(4).nV(2,3)*eye(3), ...
            Segment(4).nV(3,3)*eye(3)];
        NV44 = [Segment(4).nV(1,4)*eye(3),...
            (1 + Segment(4).nV(2,4))*eye(3), ...
            - Segment(4).nV(2,4)*eye(3), ...
            Segment(4).nV(3,4)*eye(3)];
        NV54 = [Segment(4).nV(1,5)*eye(3),...
            (1 + Segment(4).nV(2,5))*eye(3), ...
            - Segment(4).nV(2,5)*eye(3), ...
            Segment(4).nV(3,5)*eye(3)];
        NV64 = [Segment(4).nV(1,6)*eye(3),...
            (1 + Segment(4).nV(2,6))*eye(3), ...
            - Segment(4).nV(2,6)*eye(3), ...
            Segment(4).nV(3,6)*eye(3)];
        % Interpolation matrices
        Nn13 = [Segment(3).nn(1,1)*eye(3),...
            (Segment(3).nn(2,1))*eye(3), ...
            - Segment(3).nn(2,1)*eye(3), ...
            Segment(3).nn(3,1)*eye(3)];
        Nn23 = [Segment(3).nn(1,2)*eye(3),...
            (Segment(3).nn(2,2))*eye(3), ...
            - Segment(3).nn(2,2)*eye(3), ...
            Segment(3).nn(3,2)*eye(3)];
        
        % Knee joint transformation matrix (from thigh to leg)
        Joint(3).T = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(4).Q)),...
            Q2Tuv_array3(Segment(3).Q));
        % Euler angles
        Joint(3).Euler = R2mobileZXY_array3(Joint(3).T(1:3,1:3,:));
        FE = permute(Joint(3).Euler(1,1,:),[3 1 2])/pi*180;% Extension(+)/Flexion(-)
        
    otherwise
        display ('No appropriate model for the knee joint')
        
end

% -------------------------------------------------------------------------
switch Model % Complements on knee ligaments
    % Gasparutto et al. J Biomech 2015, 48: 1141–1146
    
    case {'NPN','NPS','UnPN','UnPS','UmPN','UmPS','SPN','SPS','PPN','PPS'} % Parallel knee mechanism (zero ligament length variations)
        
        % Distances (associated with segment 3)
        Joint(3).d(1,3,:) = repmat(0.0405, [1 1 n]); % ACL length
        Joint(3).d(1,4,:) = repmat(0.0433, [1 1 n]); % PCL length
        Joint(3).d(1,5,:) = repmat(0.1297, [1 1 n]); % MCL length
        
        % Initial guess for Lagrange multipliers
        lambdakK = zeros(5,1,n);
        
    case {'NLmN','NLmS','UnLmN','UnLmS','UmLmN','UmLmS','SLmN','SLmS','PLmN','PLmS'} % Minimized ligament length variations
        
        % Minimized ligament lengh variations
        Joint(3).d(1,3,:) = repmat(mean(polyval([a(1,7:-1:1),1],FE)*d3_0),[1 1 n]); % ACL
        Joint(3).d(1,4,:) = repmat(mean(polyval([a(2,7:-1:1),1],FE)*d4_0),[1 1 n]); % PCL
        Joint(3).d(1,5,:) = repmat(mean(polyval([a(3,7:-1:1),1],FE)*d5_0),[1 1 n]); % MCL
        Joint(3).d(1,6,:) = repmat(mean(polyval([a(4,7:-1:1),1],FE)*d6_0),[1 1 n]); % LCL
        
        % Weight factors on the knee ligaments kinematic constraints
        % Weight matrix
        wm = ones((size(Segment(2).rM,2) + ...
            size(Segment(3).rM,2) + ...
            size(Segment(4).rM,2) + ...
            size(Segment(5).rM,2))*3,1);
        wkK2 = [1e3;1e4;1e2;1e0]; % Weight factors on the kinematic constraints 2: ligaments
        WmkK2 = repmat(diag([wm;wkK2]),[1,1,n]);
        % The choice of weight factors can be modified
        
        % Initial guess for Lagrange multipliers
        lambdakK1 = zeros(2,1,n);
        
    case {'NLnN','NLnS','UnLnN','UnLnS','UmLnN','UmLnS','SLnN','SLnS','PLnN','PLnS'} % Targeted (nominal) ligament length variations
        
        % Targeted ligament lengh variation
        Joint(3).d(1,3,:) = permute(polyval([a(1,7:-1:1),1],FE)*d3_0,[3 2 1]); % ACL
        Joint(3).d(1,4,:) = permute(polyval([a(2,7:-1:1),1],FE)*d4_0,[3 2 1]); % PCL
        Joint(3).d(1,5,:) = permute(polyval([a(3,7:-1:1),1],FE)*d5_0,[3 2 1]); % MCL
        Joint(3).d(1,6,:) = permute(polyval([a(4,7:-1:1),1],FE)*d6_0,[3 2 1]); % LCL
        
        % Weight factors on the knee ligaments kinematic constraints
        % Weight matrix
        wm = ones((size(Segment(2).rM,2) + ...
            size(Segment(3).rM,2) + ...
            size(Segment(4).rM,2) + ...
            size(Segment(5).rM,2))*3,1);
        wkK2 = [1e4;1e4;1e4;1e4]; % Weight factors on the kinematic constraints 2: ligaments
        WmkK2 = repmat(diag([wm;wkK2]),[1,1,n]);
        % The choice of weight factors can be modified
        
        % Initial guess for Lagrange multipliers
        lambdakK1 = zeros(2,1,n);
        
end

% -------------------------------------------------------------------------
switch Model % Hip
    
    case {'NNN','UnNN','UmNN','SNN','PNN','NHnN','UnHnN','UmHnN', ...
            'SHnN','PHnN','NHmN','UnHmN','UmHmN','SHmN','PHmN', ...
            'NUnN','UnUnN','UmUnN','SUnN','PUnN','NUmN','UnUmN', ...
            'UmUmN','SUmN','PUmN','NSN','UnSN','UmSN','SSN','PSN', ...
            'NPN','UnPN','UmPN','SPN','PPN', ...
            'NLmN','UnLmN','UmLmN','SLmN','PLmN', ...
            'NLnN','UnLnN','UmLnN','SLnN','PLnN'} % No hip joint
        
        % Initial guess for Lagrange multipliers
        lambdakH = []; % To be concatenate
        
    case {'NNS','UnNS','UmNS','SNS','PNS','NHnS','UnHnS','UmHnS', ...
            'SHnS','PHnS','NHmS','UnHmS','UmHmS','SHmS','PHmS', ...
            'NUnS','UnUnS','UmUnS','SUnS','PUnS','NUmS','UnUmS', ...
            'UmUmS','SUmS','PUmS','NSS','UnSS','UmSS','SSS','PSS', ...
            'NPS','UnPS','UmPS','SPS','PPS', ...
            'NLmS','UnLmS','UmLmS','SLmS','PLmS', ...
            'NLnS','UnLnS','UmLnS','SLnS','PLnS'} % Spherical hip joint
        
        % Hip virtual marker mean coordinates (rV1 = rP4)
        % Expressed in  in (u5, rP5-rD5, w5)
        Segment(5).nV(:,1) = mean(Vnop_array3(...
            Segment(4).Q(4:6,1,:) - Segment(5).Q(4:6,1,:),...
            Segment(5).Q(1:3,1,:),...
            Segment(5).Q(4:6,1,:) - Segment(5).Q(7:9,1,:),...
            Segment(5).Q(10:12,1,:)),3);
        
        % Interpolation matrices
        NV15 = [Segment(5).nV(1,1)*eye(3),...
            (1 + Segment(5).nV(2,1))*eye(3), ...
            - Segment(5).nV(2,1)*eye(3), ...
            Segment(5).nV(3,1)*eye(3)];
        
        % Initial guess for Lagrange multipliers
        lambdakH = zeros(3,1,n);
        
    otherwise
        display ('No appropriate model for the hip joint')
        
        
end


%% ------------------------------------------------------------------------
% Run optimisation
% -------------------------------------------------------------------------

% Initial guess for Lagrange multipliers
lambdar = zeros(24,1,n); % 4 segments x 6 constraints per segment

% Initial value of the objective function
F = 1;
% Iteration number
step = 0;

% -------------------------------------------------------------------------
% Newton-Raphson
while max(permute(sqrt(sum(F.^2)),[3,2,1])) > 10e-12 && step < 20
    % The choice of tolerance and maximal number of iterations can be modified
    
    % Iteration number
    step = step + 1   % Display
    
    % Initialisation
    phik = []; % Vector of kinematic constraints
    Kk = [];  % Jacobian of kinematic constraints
    phir = []; % Vector of rigid body constraints
    Kr = []; % Jacobian of rigid body constraints
    dKlambdardQ = []; % Partial derivative of Jacobian * Lagrange multipliers
    phim = []; % Vector of motor constraints
    Km = []; % Jacobian of motor constraints
    
    % ---------------------------------------------------------------------
    switch Model % Ankle
        
        case {'NNN','NNS','NHnN','NHnS','NHmN','NHmS','NUnN','NUnS','NUmN', ...
                'NUmS','NSN','NSS','NPN','NPS',...
                'NLmN','NLnN','NLmS','NLnS'} % No ankle joint
            
            % Vector of kinematic constraints
            phikA = []; % To be concatenated
            
            % Jacobian of kinematic constraints
            KkA = [];  % To be concatenated
            % Joint structure
            Joint(2).Kk = KkA;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakAdQ = zeros(4*12,4*12,n);  % To be summed up
            
        case {'UnNN','UnNS','UnHnN','UnHnS','UnHmN','UnHmS','UnUnN', ...
                'UnUnS','UnUmN','UnUmS','UnSN','UnSS','UnPN','UnPS', ...
                'UnLmN','UnLmS','UnLnN','UnLnS',...
                'UmNN','UmNS','UmHnN','UmHnS','UmHmN','UmHmS','UmUnN', ...
                'UmUnS','UmUmN','UmUmS','UmSN','UmSS','UmPN','UmPS'...
                'UmLmN','UmLmS','UmLnN','UmLnS'} % Universal ankle joint
            
            % Vector of kinematic constraints
            % rD3 - rP2 = 0
            % w3.u2 - cos(thetaA) = 0
            phikA = [Segment(3).Q(7:9,1,:) - Segment(2).Q(4:6,1,:);...
                dot(Segment(3).Q(10:12,1,:),Segment(2).Q(1:3,1,:)) - ...
                repmat(cosd(thetaA),[1,1,n])];
            
            % Jacobian of kinematic constraints
            KkA = zeros(4,4*12,n); % Initialisation
            KkA(1:3,4:6,:) = repmat(-eye(3),[1,1,n]);
            KkA(1:3,19:21,:) = repmat(eye(3),[1,1,n]);
            KkA(4,1:3,:) = permute(Segment(3).Q(10:12,1,:),[2,1,3]); % w3'
            KkA(4,22:24,:) = permute(Segment(2).Q(1:3,1,:),[2,1,3]); % u2'
            % with transpose = permute( ,[2,1,3])
            % Joint structure
            Joint(2).Kk = KkA;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakAdQ = zeros(4*12,4*12,n); % Initialisation
            dKlambdakAdQ(1,22,:) = lambdakA(4,1,:);
            dKlambdakAdQ(2,23,:) = lambdakA(4,1,:);
            dKlambdakAdQ(3,24,:) = lambdakA(4,1,:);
            dKlambdakAdQ(22,1,:) = lambdakA(4,1,:);
            dKlambdakAdQ(23,2,:) = lambdakA(4,1,:);
            dKlambdakAdQ(24,3,:) = lambdakA(4,1,:);
            
        case {'SNN','SNS','SHnN','SHnS','SHmN','SHmS','SUnN','SUnS','SUmN', ...
                'SUmS','SSN','SSS','SPN','SPS',...
                'SLmN','SLmS','SLnN','SLnS'} % Spherical ankle joint
            
            % Vector of kinematic constraints
            % rD3 - rP2 = 0
            phikA = Segment(3).Q(7:9,1,:) - Segment(2).Q(4:6,1,:);
            
            % Jacobian of kinematic constraints
            KkA = zeros(3,4*12,n); % Initialisation
            KkA(1:3,4:6,:) = repmat(-eye(3),[1,1,n]);
            KkA(1:3,19:21,:) = repmat(eye(3),[1,1,n]);
            % Joint structure
            Joint(2).Kk = KkA;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakAdQ = zeros(4*12,4*12,n);
            
        case {'PNN','PNS','PHnN','PHnS','PHmN','PHmS','PUnN','PUnS', ...
                'PUmN','PUmS','PSN','PSS','PPN','PPS',...
                'PLmN','PLmS','PLnN','PLnS'} % Parallel ankle mechanism
            
            % Position of virtual markers
            rV12 = Mprod_array3(repmat(NV12,[1,1,n]),Segment(2).Q); % Lateral_centre_talus
            rV22 = Mprod_array3(repmat(NV22,[1,1,n]),Segment(2).Q); % TiCaL_calcaneus
            rV32 = Mprod_array3(repmat(NV32,[1,1,n]),Segment(2).Q); % CaFiL_calcaneus
            rV13 = Mprod_array3(repmat(NV13,[1,1,n]),Segment(3).Q); % Lateral_contact_tibia
            rV23 = Mprod_array3(repmat(NV23,[1,1,n]),Segment(3).Q); % TiCaL_tibia
            rV33 = Mprod_array3(repmat(NV33,[1,1,n]),Segment(3).Q); % CaFiL_fibula
            
            % Vector of kinematic constraints
            phikA = [rV13 - rV12;...
                dot((rV23 - rV22),(rV23 - rV22)) - repmat((Joint(2).d(1,1))^2,[1,1,n]);...
                dot((rV33 - rV32),(rV33 - rV32)) - repmat((Joint(2).d(1,2))^2,[1,1,n])];
            
            % Second parallel mechanism (1 spherical joint & 2 ligaments)
            KkA = zeros(5,4*12,n); % Initialisation
            KkA(1:3,1:12,:) = - repmat(NV12,[1,1,n]);
            KkA(1:3,13:24,:) = repmat(NV13,[1,1,n]);
            KkA(4,1:12,:) = - 2*Mprod_array3(permute(rV23 - rV22,[2,1,3]),repmat(NV22,[1,1,n]));
            KkA(4,13:24,:) = 2*Mprod_array3(permute(rV23 - rV22,[2,1,3]),repmat(NV23,[1,1,n]));
            KkA(5,1:12,:) = - 2*Mprod_array3(permute(rV33 - rV32,[2,1,3]),repmat(NV32,[1,1,n]));
            KkA(5,13:24,:) = 2*Mprod_array3(permute(rV33 - rV32,[2,1,3]),repmat(NV33,[1,1,n]));
            % with transpose = permute( ,[2,1,3])
            % Joint structure
            Joint(2).Kk = KkA;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakAdQ = zeros(4*12,4*12,n); % Initialisation
            dKlambdakAdQ(1:12,1:12,:) = Mprod_array3(lambdakA(4,1,:),repmat(2*NV22'*NV22,[1,1,n])) ...
                + Mprod_array3(lambdakA(5,1,:),repmat(2*NV32'*NV32,[1,1,n]));
            dKlambdakAdQ(1:12,13:24,:) = - Mprod_array3(lambdakA(4,1,:),repmat(2*NV22'*NV23,[1,1,n])) ...
                - Mprod_array3(lambdakA(5,1,:),repmat(2*NV32'*NV33,[1,1,n]));
            dKlambdakAdQ(13:24,1:12,:) = permute(dKlambdakAdQ(1:12,13:24,:),[2,1,3]); % Symetrical
            % with transpose = permute( ,[2,1,3])
            dKlambdakAdQ(13:24,13:24,:) = Mprod_array3(lambdakA(4,1,:),repmat(2*NV23'*NV23,[1,1,n])) ...
                + Mprod_array3(lambdakA(5,1,:),repmat(2*NV33'*NV33,[1,1,n]));
            
    end
    
    % ---------------------------------------------------------------------
    switch Model % Knee
        
        case {'NNN','NNS','UnNN','UnNS','UmNN','UmNS','SNN','SNS','PNN','PNS'} % No knee joint
            
            % Vector of kinematic constraints
            phikK = []; % To be concatenated
            
            % Jacobian of kinematic constraints
            KkK = []; % To be concatenated
            % Joint structure
            Joint(3).Kk = KkK;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakKdQ = zeros(4*12,4*12,n); % To be summed up
            
        case {'NHnN','NHnS','UnHnN','UnHnS','UmHnN','UmHnS','SHnN','SHnS', ...
                'PHnN','PHnS', ...
                'NHmN','NHmS','UnHmN','UnHmS','UmHmN','UmHmS','SHmN','SHmS', ...
                'PHmN','PHmS'} % Hinge knee joint
            
            % Vector of kinematic constraints
            % rD4 - rP3 = 0
            % w4.(rP3 - rD3) - L3*cos(thetaK1) = 0
            % w4.u3 - cos(thetaK2) = 0
            phikK = [Segment(4).Q(7:9,1,:) - Segment(3).Q(4:6,1,:);...
                dot(Segment(4).Q(10:12,1,:),...
                Segment(3).Q(4:6,1,:) - Segment(3).Q(7:9,1,:)) - ...
                repmat(Segment(3).L*cosd(thetaK1),[1,1,n]);...
                dot(Segment(4).Q(10:12,1,:),Segment(3).Q(1:3,1,:)) - ...
                repmat(cosd(thetaK2),[1,1,n])];
            
            % Jacobian of kinematic constraints
            % Initialisation
            KkK = zeros(5,4*12,n); % Initialisation
            KkK(1:3,16:18,:) = repmat(-eye(3),[1,1,n]);
            KkK(1:3,31:33,:) = repmat(eye(3),[1,1,n]);
            KkK(4,16:18,:) = permute(Segment(4).Q(10:12,1,:),[2,1,3]); % w4'
            KkK(4,19:21,:) = permute(-Segment(4).Q(10:12,1,:),[2,1,3]); % -w4'
            KkK(4,34:36,:) = permute(Segment(3).Q(4:6,1,:) - ...
                Segment(3).Q(7:9,1,:),[2,1,3]); % (rP3 - rD3)'
            KkK(5,13:15,:) = permute(Segment(4).Q(10:12,1,:),[2,1,3]); % w4'
            KkK(5,34:36,:) = permute(Segment(3).Q(1:3,1,:),[2,1,3]); % u3'
            % with transpose = permute( ,[2,1,3])
            % Joint structure
            Joint(3).Kk = KkK;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakKdQ = zeros(4*12,4*12,n); % Initialisation
            dKlambdakKdQ(13,34,:) = lambdakK(5,1,:);
            dKlambdakKdQ(14,35,:) = lambdakK(5,1,:);
            dKlambdakKdQ(15,36,:) = lambdakK(5,1,:);
            dKlambdakKdQ(16,34,:) = lambdakK(4,1,:);
            dKlambdakKdQ(17,35,:) = lambdakK(4,1,:);
            dKlambdakKdQ(18,36,:) = lambdakK(4,1,:);
            dKlambdakKdQ(19,34,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(20,35,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(21,36,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(34,13,:) = lambdakK(5,1,:);
            dKlambdakKdQ(35,14,:) = lambdakK(5,1,:);
            dKlambdakKdQ(36,15,:) = lambdakK(5,1,:);
            dKlambdakKdQ(34,16,:) = lambdakK(4,1,:);
            dKlambdakKdQ(35,17,:) = lambdakK(4,1,:);
            dKlambdakKdQ(36,18,:) = lambdakK(4,1,:);
            dKlambdakKdQ(34,19,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(35,20,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(36,21,:) = - lambdakK(4,1,:);
            
        case {'NUnN','NUnS','UnUnN','UnUnS','UmUnN','UmUnS','SUnN','SUnS', ...
                'PUnN','PUnS', ...
                'NUmN','NUmS','UnUmN','UnUmS','UmUmN','UmUmS','SUmN','SUmS', ...
                'PUmN','PUmS'} % Universal knee joint
            
            % Vector of kinematic constraints
            % rD4 - rP3 = 0
            % w4.(rP3 - rD3) - L3*cos(thetaK1) = 0
            phikK = [Segment(4).Q(7:9,1,:) - Segment(3).Q(4:6,1,:);...
                dot(Segment(4).Q(10:12,1,:),...
                Segment(3).Q(4:6,1,:) - Segment(3).Q(7:9,1,:)) - ...
                repmat(Segment(3).L*cosd(thetaK1),[1,1,n])];
            
            % Jacobian of kinematic constraints
            KkK = zeros(4,4*12,n); % Initialisation
            KkK(1:3,16:18,:) = repmat(-eye(3),[1,1,n]);
            KkK(1:3,31:33,:) = repmat(eye(3),[1,1,n]);
            KkK(4,16:18,:) = permute(Segment(4).Q(10:12,1,:),[2,1,3]); % w4'
            KkK(4,19:21,:) = permute(-Segment(4).Q(10:12,1,:),[2,1,3]); % -w4'
            KkK(4,34:36,:) = permute(Segment(3).Q(4:6,1,:) - ...
                Segment(3).Q(7:9,1,:),[2,1,3]); % (rP3 - rD3)'
            % with transpose = permute( ,[2,1,3])
            % Joint structure
            Joint(3).Kk = KkK;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakKdQ = zeros(4*12,4*12,n); % Initialisation
            dKlambdakKdQ(16,34,:) = lambdakK(4,1,:);
            dKlambdakKdQ(17,35,:) = lambdakK(4,1,:);
            dKlambdakKdQ(18,36,:) = lambdakK(4,1,:);
            dKlambdakKdQ(19,34,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(20,35,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(21,36,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(34,16,:) = lambdakK(4,1,:);
            dKlambdakKdQ(35,17,:) = lambdakK(4,1,:);
            dKlambdakKdQ(36,18,:) = lambdakK(4,1,:);
            dKlambdakKdQ(34,19,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(35,20,:) = - lambdakK(4,1,:);
            dKlambdakKdQ(36,21,:) = - lambdakK(4,1,:);
            
        case {'NSN','NSS','UnSN','UnSS','UmSN','UmSS','SSN','SSS','PSN','PSS'} % Spherical knee joint
            
            % Vector of kinematic constraints
            % rD4 - rP3 = 0
            phikK = Segment(4).Q(7:9,1,:) - Segment(3).Q(4:6,1,:);
            
            % Jacobian of kinematic constraints
            KkK = zeros(3,4*12,n); % Initialisation
            KkK(1:3,16:18,:) = repmat(-eye(3),[1,1,n]);
            KkK(1:3,31:33,:) = repmat(eye(3),[1,1,n]);
            % Joint structure
            Joint(3).Kk = KkK;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakKdQ = zeros(4*12,4*12,n);
            
        case {'NPN','NPS','UnPN','UnPS','UmPN','UmPS','SPN','SPS','PPN','PPS'} % Parallel knee mechanism
            
            % Position of virtual markers
            rV43 = Mprod_array3(repmat(NV43,[1,1,n]),Segment(3).Q); % Medial_contact_tibia
            rV53 = Mprod_array3(repmat(NV53,[1,1,n]),Segment(3).Q); % Lateral_contact_tibia
            rV63 = Mprod_array3(repmat(NV63,[1,1,n]),Segment(3).Q); % ACL_tibia
            rV73 = Mprod_array3(repmat(NV73,[1,1,n]),Segment(3).Q); % PCL_tibia
            rV83 = Mprod_array3(repmat(NV83,[1,1,n]),Segment(3).Q); % MCL_tibia
            rV14 = Mprod_array3(repmat(NV14,[1,1,n]),Segment(4).Q); % Medial_centre_femur
            rV24 = Mprod_array3(repmat(NV24,[1,1,n]),Segment(4).Q); % Lateral_centre_femur
            rV34 = Mprod_array3(repmat(NV34,[1,1,n]),Segment(4).Q); % ACL_femur
            rV44 = Mprod_array3(repmat(NV44,[1,1,n]),Segment(4).Q); % PCL_femur
            rV54 = Mprod_array3(repmat(NV54,[1,1,n]),Segment(4).Q); % MCL_femur
            
            % Direction of normals
            n13 = Mprod_array3(repmat(Nn13,[1,1,n]),Segment(3).Q); % Medial_tibia
            n23 = Mprod_array3(repmat(Nn23,[1,1,n]),Segment(3).Q); % Lateral_tibia
            
            % Vector of kinematic constraints
            phikK = [dot((rV14 - rV43),n13) - repmat(Joint(3).d(1,1),[1,1,n]);...
                dot((rV24 - rV53),n23) - repmat(Joint(3).d(1,2),[1,1,n]);...
                dot((rV34 - rV63),(rV34 - rV63)) - repmat((Joint(3).d(1,3))^2,[1,1,n]);...
                dot((rV44 - rV73),(rV44 - rV73)) - repmat((Joint(3).d(1,4))^2,[1,1,n]);...
                dot((rV54 - rV83),(rV54 - rV83)) - repmat((Joint(3).d(1,5))^2,[1,1,n])];
            
            % Jacobian of kinematic constraints
            KkK = zeros(5,4*12,n);
            KkK(1,13:24,:) = - Mprod_array3(permute(n13,[2,1,3]),repmat(NV43,[1,1,n])) ...
                + Mprod_array3(permute(rV14 - rV43,[2,1,3]),repmat(Nn13,[1,1,n]));
            KkK(1,25:36,:) = Mprod_array3(permute(n13,[2,1,3]),repmat(NV14,[1,1,n]));
            KkK(2,13:24,:) = - Mprod_array3(permute(n23,[2,1,3]),repmat(NV53,[1,1,n])) ...
                + Mprod_array3(permute(rV24 - rV53,[2,1,3]),repmat(Nn23,[1,1,n]));
            KkK(2,25:36,:) = Mprod_array3(permute(n23,[2,1,3]),repmat(NV24,[1,1,n]));
            KkK(3,13:24,:) = - 2*Mprod_array3(permute(rV34 - rV63,[2,1,3]),repmat(NV63,[1,1,n]));
            KkK(3,25:36,:) = 2*Mprod_array3(permute(rV34 - rV63,[2,1,3]),repmat(NV34,[1,1,n]));
            KkK(4,13:24,:) = - 2*Mprod_array3(permute(rV44 - rV73,[2,1,3]),repmat(NV73,[1,1,n]));
            KkK(4,25:36,:) = 2*Mprod_array3(permute(rV44 - rV73,[2,1,3]),repmat(NV44,[1,1,n]));
            KkK(5,13:24,:) = - 2*Mprod_array3(permute(rV54 - rV83,[2,1,3]),repmat(NV83,[1,1,n]));
            KkK(5,25:36,:) = 2*Mprod_array3(permute(rV54 - rV83,[2,1,3]),repmat(NV54,[1,1,n]));
            % with transpose = permute( ,[2,1,3])
            % Joint structure
            Joint(3).Kk = KkK;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakKdQ = zeros(4*12,4*12,n); % Initialisation
            dKlambdakKdQ(13:24,13:24,:) = - Mprod_array3(lambdakK(1,1,:),repmat(NV43'*Nn13 + Nn13'*NV43,[1,1,n])) ...
                - Mprod_array3(lambdakK(2,1,:),repmat(NV53'*Nn23 + Nn23'*NV53,[1,1,n])) ...
                + Mprod_array3(lambdakK(3,1,:),repmat(2*NV63'*NV63,[1,1,n])) ...
                + Mprod_array3(lambdakK(4,1,:),repmat(2*NV73'*NV73,[1,1,n])) ...
                + Mprod_array3(lambdakK(5,1,:),repmat(2*NV83'*NV83,[1,1,n]));
            dKlambdakKdQ(13:24,25:36,:) = Mprod_array3(lambdakK(1,1,:),repmat(Nn13'*NV14,[1,1,n])) ...
                + Mprod_array3(lambdakK(2,1,:),repmat(Nn23'*NV24,[1,1,n])) ...
                - Mprod_array3(lambdakK(3,1,:),repmat(2*NV63'*NV34,[1,1,n])) ...
                - Mprod_array3(lambdakK(4,1,:),repmat(2*NV73'*NV44,[1,1,n])) ...
                - Mprod_array3(lambdakK(5,1,:),repmat(2*NV83'*NV54,[1,1,n]));
            dKlambdakKdQ(25:36,13:24,:) = permute(dKlambdakKdQ(13:24,37:48,:),[2,1,3]); % Symetrical
            % with transpose = permute( ,[2,1,3])
            dKlambdakKdQ(25:36,25:36,:) = Mprod_array3(lambdakK(3,1,:),repmat(2*NV34'*NV34,[1,1,n])) ...
                + Mprod_array3(lambdakK(4,1,:),repmat(2*NV44'*NV44,[1,1,n])) ...
                + Mprod_array3(lambdakK(5,1,:),repmat(2*NV54'*NV54,[1,1,n]));
            
        case {'NLmN','NLmS','UnLmN','UnLmS','UmLmN','UmLmS','SLmN','SLmS','PLmN','PLmS', ...
                'NLnN','NLnS','UnLnN','UnLnS','UmLnN','UmLnS','SLnN','SLnS','PLnN','PLnS'} % Minimized and targeted (nominal) length variations
            
            % Position of virtual markers
            rV43 = Mprod_array3(repmat(NV43,[1,1,n]),Segment(3).Q); % Medial_contact_tibia
            rV53 = Mprod_array3(repmat(NV53,[1,1,n]),Segment(3).Q); % Lateral_contact_tibia
            rV63 = Mprod_array3(repmat(NV63,[1,1,n]),Segment(3).Q); % ACL_tibia
            rV73 = Mprod_array3(repmat(NV73,[1,1,n]),Segment(3).Q); % PCL_tibia
            rV83 = Mprod_array3(repmat(NV83,[1,1,n]),Segment(3).Q); % MCL_tibia
            rV93 = Mprod_array3(repmat(NV93,[1,1,n]),Segment(3).Q); % MCL_tibia
            rV14 = Mprod_array3(repmat(NV14,[1,1,n]),Segment(4).Q); % Medial_centre_femur
            rV24 = Mprod_array3(repmat(NV24,[1,1,n]),Segment(4).Q); % Lateral_centre_femur
            rV34 = Mprod_array3(repmat(NV34,[1,1,n]),Segment(4).Q); % ACL_femur
            rV44 = Mprod_array3(repmat(NV44,[1,1,n]),Segment(4).Q); % PCL_femur
            rV54 = Mprod_array3(repmat(NV54,[1,1,n]),Segment(4).Q); % MCL_femur
            rV64 = Mprod_array3(repmat(NV64,[1,1,n]),Segment(4).Q); % LCL_femur
            
            % Direction of normals
            n13 = Mprod_array3(repmat(Nn13,[1,1,n]),Segment(3).Q); % Medial_tibia
            n23 = Mprod_array3(repmat(Nn23,[1,1,n]),Segment(3).Q); % Lateral_tibia
            
            % Vector of kinematic constraint 1: contacts
            phikK1 = [dot((rV14 - rV43),n13) - Joint(3).d(1,1);...
                dot((rV24 - rV53),n23) - Joint(3).d(1,2)];
            
            % Vector of kinematic constraint 2: ligaments
            phikK2 = [dot((rV34 - rV63),(rV34 - rV63)) - Joint(3).d(1,3).^2;...
                dot((rV44 - rV73),(rV44 - rV73)) - Joint(3).d(1,4).^2;...
                dot((rV54 - rV83),(rV54 - rV83)) - Joint(3).d(1,5).^2;...
                dot((rV64 - rV93),(rV64 - rV93)) - Joint(3).d(1,6).^2];
            
            % Jacobian of kinematic constraints 1: contacts
            KkK1 = zeros(2,4*12,n);
            KkK1(1,13:24,:) = - Mprod_array3(permute(n13,[2,1,3]),repmat(NV43,[1,1,n])) ...
                + Mprod_array3(permute(rV14 - rV43,[2,1,3]),repmat(Nn13,[1,1,n]));
            KkK1(1,25:36,:) = Mprod_array3(permute(n13,[2,1,3]),repmat(NV14,[1,1,n]));
            KkK1(2,13:24,:) = - Mprod_array3(permute(n23,[2,1,3]),repmat(NV53,[1,1,n])) ...
                + Mprod_array3(permute(rV24 - rV53,[2,1,3]),repmat(Nn23,[1,1,n]));
            KkK1(2,25:36,:) = Mprod_array3(permute(n23,[2,1,3]),repmat(NV24,[1,1,n]));
            % Joint structure
            Joint(3).Kk = KkK1;
            
            % Jacobian of kinematic constraints 2: ligaments
            KkK2 = zeros(4,4*12,n);
            KkK2(1,13:24,:) = - 2*Mprod_array3(permute(rV34 - rV63,[2,1,3]),repmat(NV63,[1,1,n]));
            KkK2(1,25:36,:) = 2*Mprod_array3(permute(rV34 - rV63,[2,1,3]),repmat(NV34,[1,1,n]));
            KkK2(2,13:24,:) = - 2*Mprod_array3(permute(rV44 - rV73,[2,1,3]),repmat(NV73,[1,1,n]));
            KkK2(2,25:36,:) = 2*Mprod_array3(permute(rV44 - rV73,[2,1,3]),repmat(NV44,[1,1,n]));
            KkK2(3,13:24,:) = - 2*Mprod_array3(permute(rV54 - rV83,[2,1,3]),repmat(NV83,[1,1,n]));
            KkK2(3,25:36,:) = 2*Mprod_array3(permute(rV54 - rV83,[2,1,3]),repmat(NV54,[1,1,n]));
            KkK2(4,13:24,:) = - 2*Mprod_array3(permute(rV64 - rV93,[2,1,3]),repmat(NV93,[1,1,n]));
            KkK2(4,25:36,:) = 2*Mprod_array3(permute(rV64 - rV93,[2,1,3]),repmat(NV64,[1,1,n]));
            
            % Partial derivative of Jacobian * Lagrange multipliers
            % (kinematic constraints 1: contacts)
            dKlambdakK1dQ = zeros(4*12,4*12,n); % Initialisation
            dKlambdakK1dQ(13:24,13:24,:) = - Mprod_array3(lambdakK1(1,1,:),repmat(NV43'*Nn13 + Nn13'*NV43,[1,1,n])) ...
                - Mprod_array3(lambdakK1(2,1,:),repmat(NV53'*Nn23 + Nn23'*NV53,[1,1,n]));
            dKlambdakK1dQ(13:24,25:36,:) = Mprod_array3(lambdakK1(1,1,:),repmat(Nn13'*NV14,[1,1,n])) ...
                + Mprod_array3(lambdakK1(2,1,:),repmat(Nn23'*NV24,[1,1,n]));
            dKlambdakK1dQ(25:36,13:24,:) = permute(dKlambdakK1dQ(13:24,25:36,:),[2,1,3]); % Symetrical
            
    end
    
    % ---------------------------------------------------------------------
    switch Model % Hip
        
        case {'NNN','UnNN','UmNN','SNN','PNN','NHnN','UnHnN','UmHnN', ...
                'SHnN','PHnN','NHmN','UnHmN','UmHmN','SHmN','PHmN', ...
                'NUnN','UnUnN','UmUnN','SUnN','PUnN','NUmN','UnUmN', ...
                'UmUmN','SUmN','PUmN','NSN','UnSN','UmSN','SSN','PSN', ...
                'NPN','UnPN','UmPN','SPN','PPN', ...
                'NLmN','UnLmN','UmLmN','SLmN','PLmN', ...
                'NLnN','UnLnN','UmLnN','SLnN','PLnN'} % No hip joint
            
            % Vector of kinematic constraints
            phikH = []; % To be concatenated
            
            % Jacobian of kinematic constraints
            KkH = []; % To be concatenated
            % Joint structure
            Joint(4).Kk = KkH;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakHdQ = zeros(4*12,4*12,n); % To be summed up
            
        case {'NNS','UnNS','UmNS','SNS','PNS','NHnS','UnHnS','UmHnS', ...
                'SHnS','PHnS','NHmS','UnHmS','UmHmS','SHmS','PHmS', ...
                'NUnS','UnUnS','UmUnS','SUnS','PUnS','NUmS','UnUmS', ...
                'UmUmS','SUmS','PUmS','NSS','UnSS','UmSS','SSS','PSS', ...
                'NPS','UnPS','UmPS','SPS','PPS', ...
                'NLmS','UnLmS','UmLmS','SLmS','PLmS', ...
                'NLnS','UnLnS','UmLnS','SLnS','PLnS'} % Spherical hip joint
            
            % Vector of kinematic constraints
            % rV15 - rP4 = 0
            phikH = Mprod_array3(repmat(NV15,[1,1,n]),Segment(5).Q) - ...
                Segment(4).Q(4:6,1,:);
            
            % Jacobian of kinematic constraints
            KkH = zeros(3,4*12,n); % Initialisation
            KkH(1:3,28:30,:) = repmat(-eye(3),[1,1,n]);
            KkH(1:3,37:48,:) = repmat(NV15,[1,1,n]);
            % Joint structure
            Joint(4).Kk = KkH;
            
            % Partial derivative of Jacobian * Lagrange multipliers
            dKlambdakHdQ = zeros(4*12,4*12,n); % Initialisation
            
    end
    
    
    % ---------------------------------------------------------------------
    % Rigid body constraints and motor constraints
    for i = 2:5 % From i = 2 (Foot) to i = 5 (Pelvis)
        
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
        Kri = zeros(6,4*12,n); % Initialisation
        Kri(1:6,(i-2)*12+1:(i-2)*12+12,:) = permute(...
            [    2*ui,       vi,           wi,     zeros(3,1,n),zeros(3,1,n),zeros(3,1,n); ...
            zeros(3,1,n),    ui,      zeros(3,1,n),    2*vi,         wi,     zeros(3,1,n); ...
            zeros(3,1,n),   -ui,      zeros(3,1,n),   -2*vi,        -wi,     zeros(3,1,n); ...
            zeros(3,1,n),zeros(3,1,n),     ui,     zeros(3,1,n),     vi,         2*wi],[2,1,3]);
        % with transpose = permute( ,[2,1,3])
        % Segment structure
        Segment(i).Kr = Kri;
        
        % Partial derivative of Jacobian * Lagrange multipliers
        dKlambdaridQ = zeros(12,4*12,n); % Initialisation
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
        
        % Vector and Jacobian of motor constraints
        Kmi = zeros(size(Segment(i).rM,2)*3,4*12,n); % Initialisation
        phimi = []; % Initialisation
        for j = 1:size(Segment(i).rM,2)
            % Interpolation matrix
            NMij = [Segment(i).nM(1,j)*eye(3),...
                (1 + Segment(i).nM(2,j))*eye(3), ...
                - Segment(i).nM(2,j)*eye(3), ...
                Segment(i).nM(3,j)*eye(3)];
            % Vector of motor constraints
            phimi((j-1)*3+1:(j-1)*3+3,1,:) = Segment(i).rM(:,j,:) ...
                - Mprod_array3(repmat(NMij,[1,1,n]),Segment(i).Q);
            % Jacobian of motor contraints
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
    
    % ---------------------------------------------------------------------
    % Solution
    switch Model % Langrange multiplier or penalty based method
        
        case {'NLmN','NLmS','UnLmN','UnLmS','UmLmN','UmLmS','SLmN','SLmS','PLmN','PLmS', ...
                'NLnN','NLnS','UnLnN','UnLnS','UmLnN','UmLnS','SLnN','SLnS','PLnN','PLnS'} % Minimized and targeted (nominal) ligament length variations
            
            % Assembly
            phik = [phikA;phikK1;phikH];
            Kk = [KkA;KkK1;KkH];
            phimphikK2 = [phim;phikK2]; % Motor constraints and kinematic constraint 2: ligaments
            KmKkK2 = [Km;KkK2]; % Jacobian of motor constraints and kinematic constraint 2: ligaments
            lambdak = [lambdakA;lambdakK1;lambdakH];
            dKlambdakdQ = dKlambdakAdQ + dKlambdakK1dQ + dKlambdakHdQ;
            
            % Display errors
            Mean_phik = mean(Mprod_array3(permute(phik,[2,1,3]),phik),3)
            Mean_phir = mean(Mprod_array3(permute(phir,[2,1,3]),phir),3)
            Mean_phimphikK2 = mean(Mprod_array3(permute(phimphikK2,[2,1,3]),phimphikK2),3)
            
            % Compute dX
            % dX = inv(-dF/dx)*F(X)
            % F(X) = [KmKkK2'*W*phimphikK2m + [Kk;Kr]'*[lambdak;lambdar];[phik;phir]]
            % X = [Q;[lambdak;lambdar]]
            F = [Mprod_array3(permute(KmKkK2,[2,1,3]), Mprod_array3(WmkK2,phimphikK2)) + ...
                Mprod_array3(permute([Kk;Kr],[2,1,3]), [lambdak;lambdar]); ...
                [phik;phir]]; % with transpose = permute( ,[2,1,3])
            dFdX = [Mprod_array3(permute(KmKkK2,[2,1,3]),Mprod_array3(WmkK2,KmKkK2)) + ...
                dKlambdakdQ + dKlambdardQ, permute([Kk;Kr],[2,1,3]); ...
                [Kk;Kr],zeros(size([Kk;Kr],1),size([Kk;Kr],1),n)];
            dX = Mprod_array3(Minv_array3(-dFdX),F);
            
        otherwise
            
            % Assembly
            phik = [phikA;phikK;phikH];
            Kk = [KkA;KkK;KkH];
            lambdak = [lambdakA;lambdakK;lambdakH];
            dKlambdakdQ = dKlambdakAdQ + dKlambdakKdQ + dKlambdakHdQ;
            
            % Display errors
            Mean_phik = mean(Mprod_array3(permute(phik,[2,1,3]),phik),3)
            Mean_phir = mean(Mprod_array3(permute(phir,[2,1,3]),phir),3)
            Mean_phim = mean(Mprod_array3(permute(phim,[2,1,3]),phim),3)
            
            % Compute dX
            % dX = inv(-dF/dx)*F(X)
            % F(X) = [Km'*phim + [Kk;Kr]'*[lambdak;lambdar];[phik;phir]]
            % X = [Q;[lambdak;lambdar]]
            F = [Mprod_array3(permute(Km,[2,1,3]),phim) + ...
                Mprod_array3(permute([Kk;Kr],[2,1,3]), [lambdak;lambdar]); ...
                [phik;phir]]; % with transpose = permute( ,[2,1,3])
            dFdX = [Mprod_array3(permute(Km,[2,1,3]),Km) + ...
                dKlambdakdQ + dKlambdardQ, permute([Kk;Kr],[2,1,3]); ...
                [Kk;Kr],zeros(size([Kk;Kr],1),size([Kk;Kr],1),n)];
            dX = Mprod_array3(Minv_array3(-dFdX),F);
            
    end
    
    % ---------------------------------------------------------------------
    % Extraction from X
    Segment(2).Q = Segment(2).Q + dX(1:12,1,:);
    Segment(3).Q = Segment(3).Q + dX(13:24,1,:);
    Segment(4).Q = Segment(4).Q + dX(25:36,1,:);
    Segment(5).Q = Segment(5).Q + dX(37:48,1,:);
    
    switch Model
        case {'NNN'}
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'UnNN','UmNN','SNN','PNN'}
            lambdakA = lambdakA + dX(49:49 + size(lambdakA,1) - 1,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'NHnN','NHmN','NUnN','NUmN','NSN','NPN'}
            lambdakK = lambdakK + dX(49:49 + size(lambdakK,1) - 1,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'NLmN','NLnN'}
            lambdakK1 = lambdakK1 + dX(49:49 + size(lambdakK1,1) - 1,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'NNS'}
            lambdakH = lambdakH + dX(end - 4*6 - 2:end - 4*6,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'UnHnN','UnHmN','UnUnN','UnUmN','UnSN','UnPN', ...
                'UmHnN','UmHmN','UmUnN','UmUmN','UmSN','UmPN', ...
                'SHnN','SHmN','SUnN','SUmN','SSN','SPN', ...
                'PHnN','PHmN','PUnN','PUmN','PSN','PPN'}
            lambdakA = lambdakA + dX(49:49 + size(lambdakA,1) - 1,1,:);
            lambdakK = lambdakK + dX(49 + size(lambdakA,1):end - 4*6,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'UnLmN','UnLnN','UmLmN','UmLnN','SLmN','SLnN','PLmN','PLnN'}
            lambdakA = lambdakA + dX(49:49 + size(lambdakA,1) - 1,1,:);
            lambdakK1 = lambdakK1 + dX(49 + size(lambdakA,1):end - 4*6,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'UnNS','UmNS','SNS','PNS'}
            lambdakA = lambdakA + dX(49:49 + size(lambdakA,1) - 1,1,:);
            lambdakH = lambdakH + dX(end - 4*6 - 2:end - 4*6,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'NHnS','NHmS','NUnS','NUmS','NSS','NPS'}
            lambdakK = lambdakK + dX(49:49 + size(lambdakK,1) - 1,1,:);
            lambdakH = lambdakH + dX(end - 4*6 - 2:end - 4*6,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'NLmS','NLnS'}
            lambdakK1 = lambdakK1 + dX(49:49 + size(lambdakK1,1) - 1,1,:);
            lambdakH = lambdakH + dX(end - 4*6 - 2:end - 4*6,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        case {'UnLmS','UnLnS','UmLmS','UmLnS','SLmS','SLnS','PLmS','PLnS'}
            lambdakA = lambdakA + dX(49:49 + size(lambdakA,1) - 1,1,:);
            lambdakK1 = lambdakK1 + dX(49 + size(lambdakA,1):end - 4*6 - 3,1,:);
            lambdakH = lambdakH + dX(end - 4*6 - 2:end - 4*6,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
        otherwise
            lambdakA = lambdakA + dX(49:49 + size(lambdakA,1) - 1,1,:);
            lambdakK = lambdakK + dX(49 + size(lambdakA,1):end - 4*6 - 3,1,:);
            lambdakH = lambdakH + dX(end - 4*6 - 2:end - 4*6,1,:);
            lambdar = lambdar + dX(end - 4*6 + 1:end,1,:);
    end
    
end
