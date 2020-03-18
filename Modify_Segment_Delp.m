%  FUNCTION
% Modify_Segment_Delp.m
%__________________________________________________________________________
%
% PURPOSE
% Insert patella as segment 4 
%
% SYNOPSIS
% Segment = Modify_Segment_Delp(Segment)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
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
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas, Florent Moissent, Edouard Jouan
% September 2012
%
% Modified by Raphael Dumas
% June 2019
% Correction (definition of Segment(4).Q)
% rP4 = Mprod_array3(repmat(NV65,[1,1,n]),Segment(5).Q) - rVs14(1,1)*u4
% (not + rVs14(1,1)*u4)
%
% Modified by Raphael Dumas
% September 2019
% Parameter from Delp's model
% SCoRE applied to the modelled movement of patella w.r.t. femur
% Renamed _Delp
% Possibility of informed segment length
%
% Modified by Florent Moissenet
% March 2020
% Scaling set by Segment(i).rScale
%__________________________________________________________________________

function Segment = Modify_Segment_Delp(Segment,Model)

% Number of frames
n = size(Segment(2).rM,3);

% -------------------------------------------------------------------------
% Insert patella as segment 4 
% -------------------------------------------------------------------------

% Pelvis
Segment(6).Q = Segment(5).Q;
Segment(6).rM = Segment(5).rM;
% Thigh
Segment(5).Q = Segment(4).Q;
Segment(5).m = Segment(4).m;
Segment(5).rCs = Segment(4).rCs;
Segment(5).Is = Segment(4).Is;
Segment(5).rM = Segment(4).rM;
Segment(5).rScale = Segment(4).rScale;
Segment(5).L = Segment(4).L;
% Patella
Segment(4).Q = []; % To be defined
Segment(4).m = 0.0975; % According to the OpenSim model
Segment(4).rCs = [0;0;0]; % rCs = rP4 (origin of the SCS)
Segment(4).Is = zeros(3,3); % Patella is considered as a ponctual mass     
Segment(4).rM = []; % No associated marker

% Segment parameter of femur
if (isfield(Segment,'L') && ~isempty(Segment(5).L)) % Informed segment length
    L(5) = Segment(5).L;
else
    L(5) = mean(sqrt(sum((Segment(5).Q(4:6,1,:) - ...
        Segment(5).Q(7:9,1,:)).^2)),3); % Mean segment length
end
alpha(5) = mean(acosd(dot(Segment(5).Q(4:6,1,:) - ...
            Segment(5).Q(7:9,1,:), Segment(5).Q(10:12,1,:))./...
            sqrt(sum((Segment(5).Q(4:6,1,:) - ...
            Segment(5).Q(7:9,1,:)).^2))),3);
beta(5) = mean(acosd(dot(Segment(5).Q(10:12,1,:), ...
            Segment(5).Q(1:3,1,:))),3);
gamma(5) = mean(acosd(dot(Segment(5).Q(1:3,1,:), ...
            Segment(5).Q(4:6,1,:) - Segment(5).Q(7:9,1,:))./...
            sqrt(sum((Segment(5).Q(4:6,1,:) - ...
            Segment(5).Q(7:9,1,:)).^2))),3);
B5 = [1, L(5)*cosd(gamma(5)), cosd(beta(5)); ...
        0, L(5)*sind(gamma(5)), (cosd(alpha(5)) - cosd(beta(5))*cosd(gamma(5)))/sind(gamma(5)); ...
        0, 0, sqrt(1 - cosd(beta(5))^2 - ((cosd(alpha(5)) - cosd(beta(5))*cosd(gamma(5)))/sind(gamma(5)))^2)];
invB5 = inv(B5);

% Interpolation matrices
% Orientation of the Z axis of femur SCS in the NSCS
nZ5 = invB5(:,3); 
NVZ5 =[nZ5(1,1)*eye(3),...
       (nZ5(2,1))*eye(3), ...
       - nZ5(2,1)*eye(3), ...
       nZ5(3,1)*eye(3)]; 

% -------------------------------------------------------------------------
% Model parameters to define segment 4
% -------------------------------------------------------------------------

% Generic model
rVs65 = ([0.0026; -0.0012; 0.0000] + ... % SCoRE results for hinge in femur
    Model.geometry.T_tibia_mepic); % From Delp mechanism to mid-epicondyle (origin at rD5)
ns15 = [0;0;1]; % Axis orientation
rVs14 = [-0.0377; 0.0318; -0.0000]; % SCoRE results for hinge in patella (origin at rD4)
% rVs24 = rD4;
ns14 = ns15; % Axes of patella and femur assumed aligned in neutral posture
rVs93 = (Model.geometry.L_tibialTuber + ... Patellar tendon insertion at TT (origin at proximal endpoint)
        Model.geometry.T_tibia_mepic); % From Delp mechanism to mid-epicondyle

% Point of the hinge axis in femur (defined based on the generic model)
% Scaled by thigh (5) scaling ratio
Segment(5).nV(:,6) = [0;-1;0] + invB5*rVs65*Segment(5).rScale; % Virtual marker 6
NV65 = [Segment(5).nV(1,6)*eye(3),... Interpolation matrix
        (1 + Segment(5).nV(2,6))*eye(3),...
        - Segment(5).nV(2,6)*eye(3), ...
        Segment(5).nV(3,6)*eye(3)];
    
% Patella length (the generic model patella length is used here as the
% subject patella length is unknown: not used for scaling, only to define rP4)
L(4) = norm(Model.geometry.T_patel_patel); 
    
% -------------------------------------------------------------------------
% Definition of Segment(4).Q
% -------------------------------------------------------------------------
u4 = Vnorm_array3((Segment(3).Q(1:3,1,:) + Segment(5).Q(1:3,1,:))/2); % Mean of thigh and shank u axes
v4 = Vnorm_array3(cross(Mprod_array3(repmat(NVZ5,[1,1,n]),Segment(5).Q),u4)); % Z axis of femur SCS in the NSCS
w4 = Vnorm_array3(cross(u4,v4));
rD4 = Mprod_array3(repmat(NV65,[1,1,n]),Segment(5).Q)...
    - rVs14(1,1)*u4 ... 
    - rVs14(2,1)*v4; % Distance from the point of the hinge axis to the origin of patella SCS
rP4 = rD4 + L(4)*v4;
Segment(4).Q = [u4;rP4;rD4;w4];
