%  FUNCTION
% Modify_Segment_PS.m
%__________________________________________________________________________
%
% PURPOSE
% Insert patella as segment 4 
%
% SYNOPSIS
% Segment = Modify_Segment_PS(Segment)
%
% INPUT
% Segment (cf. data structure in user guide)
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
% Data from CT-scan and CAD 
% Renamed _GC5
%
% Modified by Raphael Dumas
% Januaray 2020
% Renamed _PS
%__________________________________________________________________________

function Segment = Modify_Segment_PS(Segment)

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
% Patella
Segment(4).Q = []; % To be defined
Segment(4).m = 0.0975; % According to the OpenSim model
Segment(4).rCs = [0;0;0]; % rCs = rP4 (origin of the SCS)
Segment(4).Is = zeros(3,3); % Patella is considered as a ponctual mass     
Segment(4).rM = []; % No associated marker

% Segment parameter of femur
L(5) = 0.4256; % From PS CT-scan
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
% Data from CT-scan and CAD of PS
% Axes of tibia, patella and femur assumed aligned in neutral posture
% -------------------------------------------------------------------------

% Virtual palpation in CT-scan/CAD model (femur)
Tf_implant2bone = [0.9999    0.0048   -0.0096   -0.0027; ...
   -0.0049    0.9999   -0.0135   -0.0043; ...
    0.0095    0.0135    0.9999   -0.0027; ...
         0         0         0    1.0000]; 

% Hinge in femur (origin at distal endpoint)
rVs65 = [0; - L(5); 0; 0] + ... Origin was at rD5
    inv(Tf_implant2bone)*...
    [([-4.585; 4.718; 24.008*(-1)] + ...  % Medial sphere of femur throchlea (left knee)
    [-5.52; 3.136; -24.34*(-1)])... % Lateral sphere of femur throchlea (left knee)
    /2000;1]; % in m
rVs65(4,:) = [];

% Hinge in patella
rVs14 = [0; L(5); 0] + rVs65 - ...
    inv(Tf_implant2bone(1:3,1:3))*[42.1; -4.2; -2.72*(-1)]/1000; % Origin of patella in femur (left knee)


% Expressed in the patella SCS
L(4) = norm([46.51; -26.48; -3.3*(-1)]/1000 - ... % Patellar tendon origin in femur (left knee)
    [42.1; -4.2; -2.72*(-1)]/1000); % Origin of patella in femur (left knee)


% Point of the hinge axis in femur
Segment(5).nV(:,6) = invB5*rVs65; % Virtual marker 6
NV65 = [Segment(5).nV(1,6)*eye(3),... Interpolation matrix
        (1 + Segment(5).nV(2,6))*eye(3),...
        - Segment(5).nV(2,6)*eye(3), ...
        Segment(5).nV(3,6)*eye(3)];
    
    
% -------------------------------------------------------------------------
% Definition of Segment(4).Q
% -------------------------------------------------------------------------

u4 = Vnorm_array3((Segment(3).Q(1:3,1,:) + Segment(5).Q(1:3,1,:))/2); % Mean of thigh and shank u axes
v4 = Vnorm_array3(cross(Mprod_array3(repmat(NVZ5,[1,1,n]),Segment(5).Q),u4)); % Z axis of femur SCS in the NSCS
w4 = Vnorm_array3(cross(u4,v4));
rP4 = Mprod_array3(repmat(NV65,[1,1,n]),Segment(5).Q) - rVs14(1,1)*u4; % Distance from the point of the hinge axis to the origin of patella SCS 
rD4 = rP4 - L(4)*v4;
Segment(4).Q = [u4;rP4;rD4;w4];
