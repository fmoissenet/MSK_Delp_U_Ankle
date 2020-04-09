% MAIN PROGRAM
% Main_Set_Geometry_SC.m
%__________________________________________________________________________
%
% PURPOSE
% Defintion of generic model data (markers, lengths, angles...)
%
% SYNOPSIS
% N/A (i.e., main program)
%
% DESCRIPTION
%
% REFERENCE
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Mprod_array3.m
% Tinv_array3.m
% Q2Twu_array3.m
% Q2Tuv_array3.m
% R2mobileZXY_array3.m
% 
% MATLAB VERSION
% Matlab R2020a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphael Dumas
% March 2020
%__________________________________________________________________________

% Number of frames
n = Model.Informed.n;


%% ------------------------------------------------------------------------
% Generic translations form reference origins and points in Delp
% All coordinates in "ground" (origin at mid-ASIS)
% Transformations between segments are only translations
% -------------------------------------------------------------------------

Femur_Origin = [-0.0707; -0.0661; 0.0835]; % Right hip joint centre (in delpjnt.jnt)
Tibia_Origin = Femur_Origin + ... 
    [-0.00525; -0.3960; 0]; % Values of knee tx and ty at 0 degree of flexion
Talus_Origin = Tibia_Origin + ...  
    [0; -0.426; 0]; % From tibia to talus origins
Calcaneous_Origin = Talus_Origin + ...
    [-0.04877; -0.04195; 0.00792]; % From talus to calcaneus origins
Toe_Origin = Calcaneous_Origin + ...
    [0.1768; -0.002; 0.00108]; % From calcaneus to toe origins

SACR = [-0.17; 0.0; 0.0]; % Directly in ground
RASIS = [0.0; 0.0; 0.128]; % Directly in ground (ASIS is on the bone!)
RPSIS = [-0.148; 0.0; 0.047]; % Directly in ground
LPSIS = [-0.148; 0.0; -0.047]; % Directly in ground
LASIS = [0.0; 0.0; -0.128]; % Directly in ground (ASIS is on the bone!)
LJC = [-0.1055; 0.072; 0]; % Directly in ground (virtual palpation on sacrum.asc)
HJC = Femur_Origin;
TROC = Femur_Origin + [-0.005; -0.042; 0.093];
LFE = Femur_Origin + [0.0; -0.404; 0.05];
MFE = Femur_Origin + [0.0; -0.404; -0.05];
LM = Tibia_Origin + [-0.005; -0.41; 0.053];
MM = Tibia_Origin + [0.006; -0.3888; -0.038];
FH = Tibia_Origin + [0.0256; -0.04157; 0.05976];
TT = Tibia_Origin + [0.039; -0.082; 0.000]; % (Delp's thesis) TT is in the bone!
% "tibial tuberosity of the bone model is too large"
TOE = Calcaneous_Origin + ...
    [-0.02; 0.018; -0.01] + ... % From calcaneus to "foot marker" origins
    [0.190; 0.0; 0.0];
B1 = Calcaneous_Origin + ...
    [-0.02; 0.018; -0.01] + ... % From calcaneus to "foot marker" origins
    [0.13; 0.0; -0.03]; % Base of metatarsal 1
B5 = Calcaneous_Origin + ...
    [-0.02; 0.018; -0.01] + ... % From calcaneus to "foot marker" origins
    [0.13; 0.0; 0.05]; % Base of metatarsal 5
HEEL = Calcaneous_Origin + ...
    [-0.02; 0.018; -0.01] + ... % From calcaneus to "foot marker" origins
    [0; 0; 0];
M1 = [0.05; -0.93; 0.045]; % Directly in ground (virtual palpation on foot.asc)
M5 = [0.007; -0.93; 0.143]; % Directly in ground (virtual palpation on foot.asc)


% Translation from lumbar joint centre to pelvis origin
t_lumbar2pelvis = [0;0;0] - LJC; % P6 is lumbar joint centre

% Translation from mid-epicondyles (enforced in sagittal plane) to tibia origin
t_midepicondyles2tibia = [Tibia_Origin(1:2,1);0] - ...
    [(LFE(1:2,1) + MFE(1:2,1))/2; 0] ;% rP3 is mid-epicondyles

% Translation from patella proximal endpoint to origin
t_proximal2patella = [0; -(0.0264 + 0.025/2); 0]; % CoM assumed at mid-endpoints
% Patella CoM + half patella heigth on Y-axis

% Translation from mid-epicondyles (enforced in sagittal plane) to talus origin
t_midepicondyles2talus = [Talus_Origin(1:2,1);0] - ...
    [(LFE(1:2,1) + MFE(1:2,1))/2; 0]; % Virtual marker rV13 is talus origin

% Translation from mid-malleoli (enforced in sagittal plane) to talus origin
t_midmalleoli2talus = [Talus_Origin(1:2,1);0] - ...
    [(LM(1:2) + MM(1:2))/2; 0]; % Virtual marker rV12 is talus origin

% Translation from mid-malleoli (enforced in sagittal plane) to calcaneous
t_midmalleoli2calcaneus = [Calcaneous_Origin(1:2,1);0] - ...
    [(LM(1:2,1) + MM(1:2,1))/2; 0]; % rP2 is mid-malleoli

% Translation form calcaneous to toe origins
t_calcaneous2toe = Toe_Origin - Calcaneous_Origin;


%% ------------------------------------------------------------------------
% Generic segment length and width
% Based on model markers
% -------------------------------------------------------------------------

Segment(2).Generic.W = norm(M1 - M5); % Between metatarsal heads
Segment(2).Generic.L = norm((M1(1:2,1) + M5(1:2,1))/2 - ... % From mid-malleolus to mid-metatarsal heads
    (LM(1:2,1) + MM(1:2,1))/2); % In sagittal plane

% Further defined subject-specific
% Segment(3).Generic.W = norm(MM -LM); % Between malleoli
% Segment(3).Generic.L = norm((LM(1:2,1) + MM(1:2,1))/2 - ... % From mid-epicondyles to mid-malleolus
%     (LFE(1:2,1) + MFE(1:2,1))/2); % In sagittal plane
% % Segment(4).Generic.W not defined
% Segment(4).Generic.L = norm(t_proximal2patella);
% Segment(5).Generic.W = norm(MFE - LFE); % Between epicondyles
% Segment(5).Generic.L = norm([(LFE(1:2,1) + MFE(1:2,1))/2; 0] - ... % From hip joint centre to mid-epicondyles 
%     HJC); % In sagittal plane

Segment(6).Generic.W = norm(LASIS - RASIS); % Between anterior iliac spines
Segment(6).Generic.L = norm([Femur_Origin(1:2,1);0] - LJC); % From LJC to midHJC

% Other information
Segment(2).Generic.alpha = 90;
Segment(2).Generic.beta = 90;
Segment(2).Generic.gamma = ...
    acosd(dot(((M1 + M5)/2 - HEEL)/norm((M1 + M5)/2 - HEEL),...
    ((LM + MM)/2 - (M1 + M5)/2)/Segment(2).Generic.L));
Segment(3).Generic.alpha = 90;
Segment(3).Generic.beta = 90;
Segment(3).Generic.gamma = 90;
Segment(4).Generic.alpha = 90;
Segment(4).Generic.beta = 90;
Segment(4).Generic.gamma = 90;
Segment(5).Generic.alpha = 90;
Segment(5).Generic.beta = 90;
Segment(5).Generic.gamma = 90;
Segment(6).Generic.alpha = 90;
Segment(6).Generic.beta = 90;
Segment(6).Generic.gamma = ...
    acosd(dot(((LASIS + RASIS)/2 - SACR)/norm((LASIS + RASIS)/2 - SACR),...
    (LJC - [Femur_Origin(1:2,1);0])/Segment(6).Generic.L));
%
Model.Generic.Height = 1.80; % A subject that is about 1.8 m

% Further defined subjec-specific
% Joint(4).Generic.d(1,1) = 55;


%% ------------------------------------------------------------------------
% Generic virtual markers
%
% The occurance of markers can differ between Segment(i).Generic.rMs and
% Segment(i).Informed.rM
% -------------------------------------------------------------------------

Segment(2).Generic.T = eye(4);
Segment(2).Generic.T(1:3,4) = (LM + MM)/2;
Segment(2).Generic.rMs(:,1) = inv(Segment(2).Generic.T)*[TOE;1];
Segment(2).Generic.rMs(:,2) = inv(Segment(2).Generic.T)*[B1;1];
Segment(2).Generic.rMs(:,3) = inv(Segment(2).Generic.T)*[B5;1];
Segment(2).Generic.rMs(:,4) = inv(Segment(2).Generic.T)*[HEEL;1];
Segment(2).Generic.rMs(:,5) = inv(Segment(2).Generic.T)*[M1;1];
Segment(2).Generic.rMs(:,6) = inv(Segment(2).Generic.T)*[M5;1];
Segment(2).Generic.rMs(4,:) = [];

Segment(3).Generic.T = eye(4);
Segment(3).Generic.T(1:3,4) = (LFE + MFE)/2;
Segment(3).Generic.rMs(:,1) = inv(Segment(3).Generic.T)*[LM;1];
Segment(3).Generic.rMs(:,2) = inv(Segment(3).Generic.T)*[MM;1];
Segment(3).Generic.rMs(:,3) = inv(Segment(3).Generic.T)*[FH;1];
% Not TT which  is on the bone !
Segment(3).Generic.rMs(4,:) = [];

Segment(5).Generic.T = eye(4);
Segment(5).Generic.T(1:3,4) = HJC;
Segment(5).Generic.rMs(:,1) = inv(Segment(5).Generic.T)*[TROC;1]; 
Segment(5).Generic.rMs(:,2) = inv(Segment(5).Generic.T)*[LFE;1];
Segment(5).Generic.rMs(:,3) = inv(Segment(5).Generic.T)*[MFE;1];
Segment(5).Generic.rMs(4,:) = [];

Segment(6).Generic.T = eye(4);
Segment(6).Generic.T(1:3,4) = LJC;
Segment(6).Generic.rMs(:,1) = inv(Segment(6).Generic.T)*[SACR;1];
Segment(6).Generic.rMs(:,2) = inv(Segment(6).Generic.T)*[RASIS;1];
Segment(6).Generic.rMs(:,3) = inv(Segment(6).Generic.T)*[RPSIS;1];
Segment(6).Generic.rMs(:,4) = inv(Segment(6).Generic.T)*[LPSIS;1];
Segment(6).Generic.rMs(:,5) = inv(Segment(6).Generic.T)*[LASIS;1];
Segment(6).Generic.rMs(4,:) = [];


%% ------------------------------------------------------------------------
% Generic ankle parameters
% -------------------------------------------------------------------------

% Subtalar axis position in Foot segment
Segment(2).Generic.rVs(:,1) = t_midmalleoli2talus;

% Subtalar axis direction in Foot segment
Segment(2).Generic.ns(:,1) = [0.781;0.600;-0.120]/norm([0.781;0.600;-0.120]);

% Ankle axis direction in Shank segment
Segment(3).Generic.ns(:,1) = [-0.105; -0.174; 0.979]/norm([-0.105; -0.174; 0.979]);

% Ankle axis position in Shank segment
% Ankle and subtalar axes superimposed (as well as Y-axes of foot and shank in neutral position)
Segment(3).Generic.rVs(:,1) = t_midepicondyles2talus;


%% ------------------------------------------------------------------------
% Subject-specific tibiofemoral parameters
% Based on informed tibiofemoral angle
% -------------------------------------------------------------------------

% Polynomial
load p_SC;

% Virtual palpation in CT-scan/CAD model (tibia)
Segment(3).Generic.T_implant2bone = [...
    0.9986    0.0143    0.0516   -0.0017; ...
   -0.0160    0.9993    0.0341   -0.0075; ...
   -0.0510   -0.0349    0.9981   -0.0021; ...
         0         0         0    1.0000];
% Virtual palpation in CT/CAD model (femur)
Segment(5).Generic.T_implant2bone = [...
    0.9980    0.0379    0.0515   -0.0017; ...
   -0.0393    0.9988    0.0274   -0.0075; ...
   -0.0504   -0.0294    0.9983   -0.0021; ...
         0         0         0    1.0000]; % Origin at distal endpoint
Segment(5).Generic.T_implant2bone = [eye(3), ...
    Segment(5).Generic.T_implant2bone(1:3,1:3)*[0; 0.4136;0]; ... % From proximal to distal endpoint in implant coordinate system
    [0, 0, 0, 1]]*Segment(5).Generic.T_implant2bone;

% Knee joint transformation matrix (from Thigh to Shank)
Joint(3).Informed.T = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(5).Informed.Q)),...
    Q2Tuv_array3(Segment(3).Informed.Q));
% Euler angles
Joint(3).Informed.Euler = R2mobileZXY_array3(Joint(3).Informed.T(1:3,1:3,:));
FE = permute(Joint(3).Informed.Euler(1,1,:),[3 1 2])/pi*180;% Extension(+)/Flexion(-)

Segment(3).Generic.rVs43 = Mprod_array3(Tinv_array3( ... % Medial_contact_tibia (in m)
    repmat(Segment(3).Generic.T_implant2bone,[1,1,n])), ... 
    [permute(polyval(pCM3x,FE),[2 3 1]); ...
    permute(polyval(pCM3y,FE),[2 3 1]); ...
    permute(polyval(pCM3z,FE),[2 3 1]); ...
    ones(1,1,n)]);
Segment(3).Generic.rVs43(4,:,:) = [];
%
Segment(3).Generic.rVs53 = Mprod_array3(Tinv_array3( ... % Lateral_contact_tibia (in m)
    repmat(Segment(3).Generic.T_implant2bone,[1,1,n])), ... 
    [permute(polyval(pCL3x,FE),[2 3 1]); ... 
    permute(polyval(pCL3y,FE),[2 3 1]); ... 
    permute(polyval(pCL3z,FE),[2 3 1]); ...
    ones(1,1,n)]);
Segment(3).Generic.rVs53(4,:,:) = [];
% 
Segment(5).Generic.rVs15 = Mprod_array3(Tinv_array3( ... % Medial_contact_femur (in m)
    repmat(Segment(5).Generic.T_implant2bone,[1,1,n])), ... 
    [permute(polyval(pCM5x,FE),[2 3 1]); ...
    permute(polyval(pCM5y,FE),[2 3 1]); ...
    permute(polyval(pCM5z,FE),[2 3 1]); ...
    ones(1,1,n)]);
Segment(5).Generic.rVs15(4,:,:) = [];
%
Segment(5).Generic.rVs25 = Mprod_array3(Tinv_array3( ... % Medial_contact_femur (in m)
    repmat(Segment(5).Generic.T_implant2bone,[1,1,n])), ... 
    [permute(polyval(pCL5x,FE),[2 3 1]); ...
    permute(polyval(pCL5y,FE),[2 3 1]); ...
    permute(polyval(pCL5z,FE),[2 3 1]); ...
    ones(1,1,n)]);
Segment(5).Generic.rVs25(4,:,:) = [];


%% ------------------------------------------------------------------------
% Subject-specific patellofemoral parameters
% -------------------------------------------------------------------------

% Hinge in femur
Segment(5).Generic.rVs(:,6) = inv(Segment(5).Generic.T_implant2bone)*...
    [([-5.0367; 4.7691; 23.377*(-1)]/1000 + ...  % Medial sphere of femur throchlea (left knee)
    [-5.937; 2.786; -23.594*(-1)]/1000)... % Lateral sphere of femur throchlea (left knee)
    /2;1]; % in m
Segment(5).Generic.rVs(4,:) = [];

% Axis orientation
Segment(5).Generic.ns(:,1) = inv(Segment(5).Generic.T_implant2bone(1:3,1:3))*...
    ([-5.937; 2.786; -23.594*(-1)] - [-5.0367; 4.7691; 23.377*(-1)])/...
    norm([-5.937; 2.786; -23.594*(-1)] - [-5.0367; 4.7691; 23.377*(-1)]);

% Hinge in patella
Segment(4).Generic.rVs(:,1) = inv(Segment(5).Generic.T_implant2bone(1:3,1:3))*...
        (([-5.0367; 4.7691; 23.377*(-1)]/1000 + ...  % Medial sphere of femur throchlea (left knee)
    [-5.937; 2.786; -23.594*(-1)]/1000)/2 -... % Lateral sphere of femur throchlea (left knee)
    [38.0; -1.53; -1.23*(1)]/1000); % Origin of patella (left knee)

% Axes of patella and femur assumed aligned in neutral posture
Segment(4).Generic.ns(:,1) = Segment(5).Generic.ns(:,1); 

% Patellar tendon origin
% rVs24 = rD4 
% Patellar tendon insertion
Segment(3).Generic.rVs(1:4,9) = inv(Segment(3).Generic.T_implant2bone)*...
    [[34.7; -72.53; -6.14*(-1)]/1000;1]; % Patellar tendon insertion (left knee)
Segment(3).Generic.rVs(4,:) = [];

   
%% ------------------------------------------------------------------------
% Generic hip parameters
% -------------------------------------------------------------------------
    
% Hip
Segment(6).Generic.rVs(:,1) = t_lumbar2pelvis + Femur_Origin; % Hip joint centre on Z-axis


%% ------------------------------------------------------------------------
% Generic muscle parameters
% -------------------------------------------------------------------------

% Origin / insertion in Segment 2
% Virtual markers start at 4 (others are for ankle joint)
%
% Insertion of gastrocnemius medialis
Segment(2).Generic.rVs(:,4)  = t_midmalleoli2calcaneus + ...
    [0.00440000; 0.03100000; -0.00530000];
% Insertion of gastrocnemius lateralis
Segment(2).Generic.rVs(:,5)  = t_midmalleoli2calcaneus + ...
    [0.00440000; 0.03100000; -0.00530000];
% Insertion of soleus
Segment(2).Generic.rVs(:,6) = t_midmalleoli2calcaneus + ...
    [0.00440000; 0.03100000; -0.00530000];
% Via point of tibialis posterior (considered as insertion)
Segment(2).Generic.rVs(:,7)  = t_midmalleoli2calcaneus + ...
    [0.04170000; 0.03340000; -0.02860000];
% Insertion of tibialis anterior
Segment(2).Generic.rVs(:,8) = t_midmalleoli2calcaneus + ...
    [0.11660000; 0.01780000; -0.03050000];
% Via point of peroneus brevis (considered as insertion)
Segment(2).Generic.rVs(:,9)  = t_midmalleoli2calcaneus + ...
    [0.04710000; 0.02700000; 0.02330000];
% Via point of peroneus longus (considered as insertion)
Segment(2).Generic.rVs(:,10) = t_midmalleoli2calcaneus + ...
    [0.04380000; 0.02300000; 0.02210000];
% Insertion of peroneus tertius
Segment(2).Generic.rVs(:,11) = t_midmalleoli2calcaneus + ...
    [0.08570000; 0.02280000; 0.02990000];
% Via point of extensor digitorum longus (considered as insertion)
Segment(2).Generic.rVs(:,12) = t_midmalleoli2calcaneus + ...
    [0.09220000; 0.03880000; -0.00010000];
% Via point of extensor hallucis longus (considered as insertion)
Segment(2).Generic.rVs(:,13) = t_midmalleoli2calcaneus + ...
    [0.09700000; 0.03890000; -0.02110000];
% Via point of flexor digitorum longus (considered as insertion)
Segment(2).Generic.rVs(:,14) = t_midmalleoli2calcaneus + ...
    [0.04360000; 0.03150000; -0.02800000];
% Via point of flexor hallucis longus (considered as insertion)
Segment(2).Generic.rVs(:,15) = t_midmalleoli2calcaneus + ...
    [0.03740000; 0.02760000; -0.02410000];

% Origin / insertion in Segment 3
% Virtual markers start at 10
% (others are for ankle, tibio-femoral and patello-femoral joints)
%
% Insertion of tensor fasciae latae
Segment(3).Generic.rVs(:,10) = t_midepicondyles2tibia + ...
    [0.0060; -0.0487; 0.0297];
% Via point of gracillis (considered as insertion)
Segment(3).Generic.rVs(:,11) = t_midepicondyles2tibia + ...
    [-0.0154; -0.0475; -0.0358];
% Insertion of semimembranosus
Segment(3).Generic.rVs(:,12) = t_midepicondyles2tibia + ...
    [-0.0243; -0.0536; -0.0194];
% Via point of semitendinus (considered as insertion)
Segment(3).Generic.rVs(:,13) = t_midepicondyles2tibia + ...
    [-0.0314; -0.0545; -0.0146];
% Insertion of biceps femoris long head
Segment(3).Generic.rVs(:,14) = t_midepicondyles2tibia + ...
    [-0.0081; -0.0729; 0.0423];
% Via point of sartorius (considered as insertion)
Segment(3).Generic.rVs(:,15) = t_midepicondyles2tibia + ...
    [-0.0056; -0.0419; -0.0399];
% Insertion of biceps femoris short head
Segment(3).Generic.rVs(:,16) = t_midepicondyles2tibia + ...
    [-0.0101; -0.0725; 0.0406];
% Via point of gastrocnemius medialis
Segment(3).Generic.rVs(:,17) = t_midepicondyles2tibia + ...
    [-0.0217; -0.0487; -0.0295];
% Via point of gastrocnemius lateralis
Segment(3).Generic.rVs(:,18) = t_midepicondyles2tibia + ...
    [-0.0242; -0.0481; 0.0235];
% Origin of soleus
Segment(3).Generic.rVs(:,19) = t_midepicondyles2tibia + ...
    [-0.0024; -0.1533; 0.0071];
% Via point of tibialis posterior (considered as origin)
Segment(3).Generic.rVs(:,20) = t_midepicondyles2tibia + ...
    [-0.0144; -0.4051; -0.0229];
% Via point of tibialis anterior (considered as origin)
Segment(3).Generic.rVs(:,21) = t_midepicondyles2tibia + ...
    [0.0329; -0.3951; -0.0177];
% Via point of peroneus brevis (considered as origin)
Segment(3).Generic.rVs(:,22) = t_midepicondyles2tibia + ...
    [-0.0144; -0.4295; 0.0289];
% Via point of peroneus longus (considered as origin)
Segment(3).Generic.rVs(:,23) = t_midepicondyles2tibia + ...
    [-0.0162; -0.4319; 0.0289];
% Via point of peroneus tertius (considered as origin)
Segment(3).Generic.rVs(:,24) = t_midepicondyles2tibia + ...
    [0.0229; -0.4069; 0.0159];
% Via point of extensor digitorum longus (considered as origin)
Segment(3).Generic.rVs(:,25) = t_midepicondyles2tibia + ...
    [0.0289; -0.4007; 0.0072];
% Via point of extensor hallucis longus (considered as origin)
Segment(3).Generic.rVs(:,26) = t_midepicondyles2tibia + ...
    [0.0326; -0.3985; -0.0085];
% Via point of flexor digitorum longus (considered as origin)
Segment(3).Generic.rVs(:,27) = t_midepicondyles2tibia + ...
    [-0.0154; -0.4051; -0.0196];
% Via point of flexor hallucis longus (considered as origin)
Segment(3).Generic.rVs(:,28) = t_midepicondyles2tibia + ...
    [-0.0186; -0.4079; -0.0174];

% Origin / insertion in Segment 4
% Virtual markers start at 2 (other is for patello-femoral joint)
% Origin at distal endpoint
%
% Insertion of rectus femoris
Segment(4).Generic.rVs(:,2) = t_proximal2patella + ...
    [0.0121; 0.0437; -0.0010];
% Insertion of vastus medialis
Segment(4).Generic.rVs(:,3) = t_proximal2patella + ...
    [0.0063; 0.0445; -0.0170];
% Insertion of vastus intermedialis
Segment(4).Generic.rVs(:,4) = t_proximal2patella + ...
    [0.0058; 0.0480; -0.0006];
% Insertion of vastus lateralis
Segment(4).Generic.rVs(:,5) = t_proximal2patella + ...
    [0.0103; 0.0423; 0.0141];

% Origin / insertion in Segment 5
% Virtual markers start at 7
% (others are for tibiofemoral and patellofemoral joint)
%
% Via point of gluteus maximus I (considered as insertion)
Segment(5).Generic.rVs(:,7) = [-0.0457; -0.0248; 0.0392];
% Via point of gluteus maximus II (considered as insertion)
Segment(5).Generic.rVs(:,8) = [-0.0426; -0.0530; 0.0293];
% Via point of gluteus maximus III (considered as insertion)
Segment(5).Generic.rVs(:,9) = [-0.0299; -0.1041; 0.0135];
% Insertion of gluteus medius I
Segment(5).Generic.rVs(:,10) = [-0.0218; -0.0117; 0.0555];
% Insertion of gluteus medius II
Segment(5).Generic.rVs(:,11) = [-0.0258; -0.0058; 0.0527]; 
% Insertion of gluteus medius III
Segment(5).Generic.rVs(:,12) = [-0.0309; -0.0047; 0.0518];
% Insertion of gluteus minimus I
Segment(5).Generic.rVs(:,13) = [-0.0072; -0.0104; 0.0560];
% Insertion of gluteus minimus II
Segment(5).Generic.rVs(:,14) = [-0.0096; -0.0104; 0.0560];
% Insertion of gluteus minimus III
Segment(5).Generic.rVs(:,15) = [-0.0135; -0.0083; 0.0550];
% Insertion of adductor longus
Segment(5).Generic.rVs(:,16) = [0.0050; -0.2111; 0.0234];
% Insertion of adductor brevis
Segment(5).Generic.rVs(:,17) = [0.0009; -0.1196; 0.0294];
% Insertion of adductor magnus I
Segment(5).Generic.rVs(:,18) = [-0.0045; -0.1211; 0.0339];
% Insertion of adductor magnus II
Segment(5).Generic.rVs(:,19) = [0.0054; -0.2285; 0.0227];
% Insertion of adductor magnus III
Segment(5).Generic.rVs(:,20) = [0.0070; -0.3837; -0.0266];  
% Insertion of pectineus
Segment(5).Generic.rVs(:,21) = [-0.0122; -0.0822; 0.0253];
% Via point of illiacus (considered as insertion)
Segment(5).Generic.rVs(:,22) = [ 0.0017; -0.0543; 0.0057];
% Via point of psoas (considered as insertion)
Segment(5).Generic.rVs(:,23) = [ 0.0016; -0.0507; 0.0038]; 
% Insertion of quadratus femoris
Segment(5).Generic.rVs(:,24) = [-0.0381; -0.0359; 0.0366]; 
% Insertion of gemellus
Segment(5).Generic.rVs(:,25) = [-0.0142; -0.0033; 0.0443]; 
% Insertion of piriformis
Segment(5).Generic.rVs(:,26) = [-0.0148; -0.0036; 0.0437];
% Via point tensor fasciae latae (most proximal on thigh)
Segment(5).Generic.rVs(:,27) = [0.0294; -0.0995; 0.0597];
% Via point tensor fasciae latae (most distal on thigh)
Segment(5).Generic.rVs(:,28) = [0.0054; -0.4049; 0.0357];
% Via point of sartorius
Segment(5).Generic.rVs(:,29) = [-0.0030; -0.3568; -0.0421];
% Origin of biceps femoris short head
Segment(5).Generic.rVs(:,30) = [0.0050; -0.2111; 0.0234];

% % Via point of rectus femoris (from -83.65 to -150 degrees of knee flexion)
% Segment(5).Generic.rVs(:,31) = [0.0334; -0.4030; 0.0019];

% Via point of vastus medialis (considered as origin from 0 to -69.32 degrees of knee flexion)
Segment(5).Generic.rVs(:,32) = [0.0356; -0.2769; 0.0009];
% % Via point of vastus medialis (considered as origin from -69.32 to -101.99 degrees of knee flexion)
% Segment(5).Generic.rVs(:,33) = [0.0370; -0.4048; -0.0125];
% % Via point of vastus medialis (considered as origin from -101.99 to -150 of knee flexion)
% Segment(5).Generic.rVs(:,34) = [0.0274; -0.4255; -0.0131];

% Via point of vastus intermedialis (considered as origin from 0 to -81.36 degrees of knee flexion)
Segment(5).Generic.rVs(:,35) = [0.0335; -0.2084; 0.0285];
% % Via point of vastus intermedialis (considered as origin from -81.36 to -150 degrees of knee flexion)
% Segment(5).Generic.rVs(:,36) = [0.0343; -0.4030; 0.0055];

% Via point of vastus lateralis (considered as origin from 0 to -69.32 degrees of knee flexion)
Segment(5).Generic.rVs(:,37) = [0.0269; -0.2591; 0.0409];
% % Via point of vastus lateralis (considered as origin from -69.32 to -110 degrees of knee flexion)
% Segment(5).Generic.rVs(:,38) = [0.0361; -0.4030; 0.0205];
% % Via point of vastus lateralis (considered as origin from -110 to -150 degrees of knee flexion)
% Segment(5).Generic.rVs(:,39) = [0.0253; -0.4243; 0.0184];

% Via point (considered as origin) of gastrocnemius medialis (from 5.73 to -44.12 degrees of knee flexion)
Segment(5).Generic.rVs(:,40) = [-0.0239; -0.4022; -0.0258];
% Origin of gastrocnemius medialis (from -44.12 to -150 degrees of knee flexion)
Segment(5).Generic.rVs(:,41) = [-0.0127; -0.3929; -0.0235];
% Virtual points V40 and V41 are both required to compute MT lenght

% Via point (considered as origin) of gastrocnemius lateralis (from 5.73 to -44.12 degrees of knee flexion)
Segment(5).Generic.rVs(:,42) = [-0.0254; -0.4018; 0.0274];
% Origin of gastrocnemius lateralis (from -44.12 to -150 degrees of knee flexion)
Segment(5).Generic.rVs(:,43) = [-0.0155; -0.3946; 0.0272];
% Virtual points V42 and V43 are both required to compute MT lenght

% Origin / insertion in Segment 6
% Virtual markers start at 2 (other is for hip joint)
%
% Via point gluteus maximus I (considered as origin)
Segment(6).Generic.rVs(:,2) = t_lumbar2pelvis + ...
    [-0.1291; 0.0012; 0.0886];
% Via point gluteus maximus II (considered as origin)
Segment(6).Generic.rVs(:,3) = t_lumbar2pelvis + ...
    [-0.1376; -0.0520; 0.0914]; 
% Via point gluteus maximus III (considered as origin)
Segment(6).Generic.rVs(:,4) = t_lumbar2pelvis + ...
    [-0.1529;-0.1052; 0.0403];
% Origin of gluteus medius I
Segment(6).Generic.rVs(:,5) = t_lumbar2pelvis + ...
    [-0.0408; 0.0304; 0.1209];
% Origin of gluteus medius II
Segment(6).Generic.rVs(:,6) = t_lumbar2pelvis + ...
    [-0.0855; 0.0445; 0.0766];
% Origin of gluteus medius III
Segment(6).Generic.rVs(:,7) = t_lumbar2pelvis + ...
    [-0.1223; 0.0105; 0.0648];
% Origin of gluteus minimus I
Segment(6).Generic.rVs(:,8) = t_lumbar2pelvis + ...
    [-0.0467;-0.0080; 0.1056];
% Origin of gluteus minimus II
Segment(6).Generic.rVs(:,9) = t_lumbar2pelvis + ...
    [-0.0633;-0.0065; 0.0991];
% Origin of gluteus minimus III
Segment(6).Generic.rVs(:,10) = t_lumbar2pelvis + ...
    [-0.0834;-0.0063; 0.0856];
% Origin of adductor longus
Segment(6).Generic.rVs(:,11) = t_lumbar2pelvis + ...
    [-0.0316;-0.0836; 0.0169];
% Origin of adductor brevis
Segment(6).Generic.rVs(:,12) = t_lumbar2pelvis + ...
    [-0.0587;-0.0915; 0.0164];
% Origin of adductor magnus I
Segment(6).Generic.rVs(:,13) = t_lumbar2pelvis + ...
    [-0.0732;-0.1174; 0.0255];
% Origin of adductor magnus II
Segment(6).Generic.rVs(:,14) = t_lumbar2pelvis + ...
    [-0.0831;-0.1192; 0.0308];
% Origin of adductor magnus III
Segment(6).Generic.rVs(:,15) = t_lumbar2pelvis + ...
    [-0.0771;-0.1181; 0.0276];
% Origin of pectineus
Segment(6).Generic.rVs(:,16) = t_lumbar2pelvis + ...
    [-0.0431;-0.0768; 0.0451];
% Via point of illiacus (considered as origin)
Segment(6).Generic.rVs(:,17) = t_lumbar2pelvis + ...
    [-0.0218;-0.0550; 0.0851];
% Via point of psoas (considered as origin)
Segment(6).Generic.rVs(:,18) = t_lumbar2pelvis + ...
    [-0.0238;-0.0570; 0.0759];
% Origin of quadrutus femoris
Segment(6).Generic.rVs(:,19) = t_lumbar2pelvis + ...
    [-0.1143;-0.1151; 0.0520];
% Origin of gemellus
Segment(6).Generic.rVs(:,20) = t_lumbar2pelvis + ...
    [-0.1133;-0.0820; 0.0714];
% Via point piriformis
Segment(6).Generic.rVs(:,21) = t_lumbar2pelvis + ...
    [-0.1193;-0.0276; 0.0657];
% Origin of tensor fasciae latae
Segment(6).Generic.rVs(:,22) = t_lumbar2pelvis + ...
    [-0.0311; 0.0214; 0.1241];
% Origin of gracillis
Segment(6).Generic.rVs(:,23) = t_lumbar2pelvis + ...
    [-0.0563;-0.1038; 0.0079];
% Origin of semimembranosus
Segment(6).Generic.rVs(:,24) = t_lumbar2pelvis + ...
    [-0.1192;-0.1015; 0.0695];
% Origin of semitendinus
Segment(6).Generic.rVs(:,25) = t_lumbar2pelvis + ...
    [-0.1237;-0.1043; 0.0603];
% Origin of biceps femoris long head
Segment(6).Generic.rVs(:,26) = t_lumbar2pelvis + ...
    [-0.1244;-0.1001; 0.0666];
% Origin of sartorius
Segment(6).Generic.rVs(:,27) = t_lumbar2pelvis + ...
    [-0.0153;-0.0013; 0.1242];
% Origin of rectus femoris
Segment(6).Generic.rVs(:,28) = t_lumbar2pelvis + ...
    [-0.0295;-0.0311; 0.0968];


% -------------------------------------------------------------------------
% Constant lenght between via points embedded in the same segment
% -------------------------------------------------------------------------

% Gluteus maximus I
Model.Generic.Lcst(1,:) = [0, 0, 0, 0 ... % No lenght embedded in segment 1, 2, 3, or 4
        norm([-0.0457;-0.0248; 0.0392] - ... % Via point 2 of gluteus maximus I (in thigh)
    [-0.0277; -0.0566; 0.0470]), ... % Insertion of gluteus maximus I (in thigh)
    norm([-0.1195; 0.0612; 0.0700] - ... % Origin of gluteus maximus I (in pelvis)
    [-0.1291; 0.0012; 0.0886])]; % Via point 1 of gluteus maximus I (in pelvis)

% Gluteus maximus II
Model.Generic.Lcst(2,:) = [0, 0, 0, 0 ... % No lenght embedded in segment 1, 2, 3, or 4
    norm([-0.0426;-0.0530; 0.0293] - ... % Via point 2 of gluteus maximus II (in thigh)
    [-0.0156; -0.1016; 0.0419]), ... % Insertion of gluteus maximus II (in thigh)
    norm([-0.1349; 0.0176; 0.0563] - ... % Origin of gluteus maximus II (in pelvis)
    [-0.1376;-0.0520; 0.0914])]; % Via point 1 of gluteus maximus II (in pelvis)

% Gluteus maximus III
Model.Generic.Lcst(3,:) = [0, 0, 0, 0 ... % No lenght embedded in segment 1, 2, 3, or 4
    norm([-0.0299;-0.1041; 0.0135] - ... % Via point 2 of gluteus maximus III (in thigh)
    [-0.0060; -0.1419; 0.0411]), ... % Insertion of gluteus maximus III (in thigh)
    norm([-0.1556; -0.0314; 0.0058] - ... % Origin of gluteus maximus III (in pelvis)
    [-0.1529;-0.1052; 0.0403])]; % Via point 1 of gluteus maximus III (in pelvis)

% Gluteus medius I
Model.Generic.Lcst(4,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Gluteus medius II
Model.Generic.Lcst(5,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Gluteus medius III
Model.Generic.Lcst(6,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Gluteus minimus I
Model.Generic.Lcst(7,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Gluteus minimus II
Model.Generic.Lcst(8,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Gluteus minimus III
Model.Generic.Lcst(9,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Adductus longus
Model.Generic.Lcst(10,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Adductus brevis
Model.Generic.Lcst(11,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Adductus magnus I
Model.Generic.Lcst(12,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Adductus magnus II
Model.Generic.Lcst(13,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Adductus magnus III
Model.Generic.Lcst(14,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Pectineus
Model.Generic.Lcst(15,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6

% Iliacus
Model.Generic.Lcst(16,:) = [0, 0, 0, 0 ... % No lenght embedded in segment 1, 2, 3, or 4
    norm([0.0017;-0.0543; 0.0057] - ... % Via point 2 of illiacus (in thigh)
    [-0.0193; -0.0621; 0.0129]), ... % Insertion of illiacus (in thigh)
    norm([-0.0674; 0.0365; 0.0854] - ... % Origin of illiacus (in pelvis)
    [-0.0218;-0.0550; 0.0851])]; % Via point 1 of illiacus (in pelvis)

% Psoas
Model.Generic.Lcst(17,:) = [0, 0, 0, 0 ... % No lenght embedded in segment 1, 2, 3, or 4
    norm([0.0016;-0.0507; 0.0038] - ... % Via point 2 of psoas (in thigh)
    [-0.0188; -0.0597; 0.0104]), ... % Insertion of psoas (in thigh)
    norm([-0.0647; 0.0887; 0.0289] - ... % Origin of psoas (in pelvis)
    [-0.0238;-0.0570; 0.0759])]; % Via point 1 of psoas (in pelvis)

% Quadratus femoris
Model.Generic.Lcst(18,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Gemelli
Model.Generic.Lcst(19,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6

% Piriformis
Model.Generic.Lcst(20,:) = [0, 0, 0, 0, 0 ... % No lenght embedded in segment 1, 2, 3, 4, or 5
    norm([-0.1396; 0.0003; 0.0235] - ... % Origin of piriformis (in pelvis)
    [-0.1193;-0.0276; 0.0657])]; % Via point 1 of piriformis (in pelvis)

% Tensor fascia latae
Model.Generic.Lcst(21,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6

% Gracilis
Model.Generic.Lcst(22,:) = [0, 0, ... % No lenght embedded in segment 1 or 2
    norm([-0.0154;-0.0475;-0.0358] - ... % Via point 1 of gracilis (in shank)
    [0.006; -0.0836; -0.0228]), ... % Insertion of gracilis (in shank)
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Sartorius
Model.Generic.Lcst(23,:) = [0, 0, ... % No lenght embedded in segment 1 or 2
    norm([-0.0056;-0.0419;-0.0399] - ... % Via point 2 of sartorius (in shank)
    [0.006; -0.0589; -0.0383]) + ... % Via point 3 of sartorius (in shank)
    norm([0.006; -0.0589; -0.0383] - ... % Via point 3 of sartorius (in shank)
    [0.0243; -0.0840; -0.0252]), ... % Insertion of sartorius (in shank)
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Semimembranosus
Model.Generic.Lcst(24,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6

% Semitendinosus
Model.Generic.Lcst(25,:) = [0, 0, ... % No lenght embedded in segment 1 or 2
    norm([-0.0314;-0.0545;-0.0146] - ... % Via point 1 of semitendinosus (in shank)
    [-0.0113; -0.0746; -0.0245]) + ... % Via point 2 of semitendinosus (in shank)
    norm([-0.0113; -0.0746; -0.0245] - ... % Via point 2 of semitendinosus (in shank)
    [0.0027; -0.0956; -0.0193]), ... % Insertion of semitendinosus (in shank)
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Biceps femoris long head
Model.Generic.Lcst(26,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Biceps femoris short head
Model.Generic.Lcst(27,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6

% Rectus femoris
% From 0 to -83.65 and from -83.65 to -150 degrees of knee flexion
Model.Generic.Lcst(28,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6

% Vastus medialis
% From 0 to -69.32 degrees of knee flexion and beyond
Model.Generic.Lcst(29,:) = [0, 0, 0, 0, ... % No lenght embedded in segment 1, 2, 3, or 4
    norm([0.0140; -0.2099; 0.0188] - ... % Origin of vastus medialis (in thigh)
    [0.0356;-0.2769; 0.0009]), ... % Via point 1 of vastus medialis (from 0 to -69.32 degrees of knee flexion) (in thigh)
    0]; % No lenght embedded in segment 6

% Vastus intermedialis
% From 0 to -81.36 degrees of knee flexion and beyond
Model.Generic.Lcst(30,:) = [0, 0, 0, 0, ... % No lenght embedded in segment 1, 2, 3, or 4
    norm([0.0290; -0.1924; 0.0310] - ... % Origin of vastus intermedialis (in thigh)
    [0.0335;-0.2084; 0.0285]), ... % Via point 1 of vastus intermedialis (from 0 to -81.36 degrees of knee flexion) (in thigh)
    0]; % No lenght embedded in segment 6

% Vastus lateralis
% From 0 to -69.32 degrees of knee flexion and beyond
Model.Generic.Lcst(31,:) = [0, 0, 0, 0, ... % No lenght embedded in segment 1, 2, 3, or 4
    norm([0.0048; -0.1854; 0.0349] - ... % Origin of vastus lateralis (in thigh)
    [0.0269;-0.2591; 0.0409]), ... % Via point 1 of vastus lateralis (from 0 to -69.32 degrees of knee flexion) (in thigh)
    0]; % No lenght embedded in segment 6

% Gastrocnemius medialis
% From 0 to -44.14 degrees of knee flexion and beyond
Model.Generic.Lcst(32,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Gastrocnemius lateralis
% From 0 to -44.14 degrees of knee flexion and beyond
Model.Generic.Lcst(33,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6
% Soleus
Model.Generic.Lcst(34,:) = zeros(1,6); % No lenght embedded in segment 1, 2, 3, 4, 5, or 6

% Tibialis posterior
Model.Generic.Lcst(35,:) = [0, ... % No lenght embedded in segment 1
    norm([0.04170000; 0.03340000; -0.02860000] - ... % Via point 2 of tibialis posterior (in foot)
    [0.07720000; 0.01590000; -0.02810000]), ... % Insertion of tibialis posterior (in foot)
    norm([-0.0094; -0.1348; 0.0019] - ... % Origin of tibialis posterior (in shank)
    [-0.0144;-0.4051;-0.0229]), ... % Via point 1 of tibialis posterior (in shank)
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Tibialis anterior
Model.Generic.Lcst(36,:) = [0, 0, ... % No lenght embedded in segment 1 or 2
    norm([0.0179; -0.1624; 0.0115] - ... % Origin of tibialis anterior (in shank)
    [0.0329;-0.3951;-0.0177]), ... % Via point 1 of tibialis anterior (in shank)
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Peroneus brevis
Model.Generic.Lcst(37,:) = [0, ... % No lenght embedded in segment 1
       norm([0.04710000; 0.02700000; 0.02330000] - ... % Via point 4 of peroneus brevis (in foot)
    [0.06770000; 0.02190000; 0.03430000]), ... % Insertion of peroneus brevis (in foot)
    norm([-0.0070; -0.2646; 0.0325] - ... % Origin of peroneus brevis (in shank)
    [-0.0198; -0.4184; 0.0283]) +  ... % Via point 1 of peroneus brevis (in shank)
    norm([-0.0198; -0.4184; 0.0283] -  ... % Via point 1 of peroneus brevis (in shank)
    [-0.0144;-0.4295; 0.0289]), ... % Via point 2 of peroneus brevis (in shank)
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Peroneus longus
Model.Generic.Lcst(38,:) = [0, ... % No lenght embedded in segment 1
     norm([0.04380000; 0.02300000; 0.02210000] - ... % Via point 3 of peroneus longus (in foot)
    [0.06810000; 0.01060000; 0.02840000]) + ... % Via point 4 of peroneus longus (in foot)
    norm([0.06810000; 0.01060000; 0.02840000] - ... % Via point 4 of peroneus longus (in foot)
    [0.08520000; 0.00690000; 0.01180000]) + ... % Via point 5 of peroneus longus (in foot)
    norm([0.08520000; 0.00690000; 0.01180000] - ... % Via point 5 of peroneus longus (in foot)
    [0.12030000; 0.00850000; -0.01840000]), ... % Insertion of peroneus longus (in foot)
    norm([0.0005; -0.1568; 0.0362] - ... % Origin of peroneus longus (in shank)
    [-0.0207; -0.4205; 0.0286]) + ... % Via point 1 of peroneus longus (in shank)
    norm([-0.0207; -0.4205; 0.0286] - ... % Via point 1 of peroneus longus (in shank)
    [-0.0162;-0.4319; 0.0289]), ... % Via point 2 of peroneus longus (in shank)
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Peroneus tertius
Model.Generic.Lcst(39,:) = [0, 0, ... % No lenght embedded in segment 1 or 2
    norm([0.0010; -0.2804; 0.0231] - ... % Origin of peroneus tertius
    [0.0229;-0.4069; 0.0159]), ... % Via point 1 of peroneus tertius
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Extensor digitorum longus
Model.Generic.Lcst(40,:) = [0, ... % No lenght embedded in segment 1
    norm([0.09220000; 0.03880000; -0.00010000] - ... % Via point 2 of extensor digitorum longus
    [0.16160000; 0.00550000; 0.01300000]) + ... % Via point 3 of extensor digitorum longus (in foot: calcaneous)
    norm([0.16160000; 0.00550000; 0.01300000] - ... % Via point 3 of extensor digitorum longus (in foot: calcaneous)
    (t_calcaneous2toe + [0.00030000; 0.00470000; 0.01530000])) + ... % Via point 4 of extensor digitorum longus (in foot: toe)
    norm([0.00030000; 0.00470000; 0.01530000] - ... % Via point 4 of extensor digitorum longus (in foot: toe)
    [0.04430000; -0.00040000; 0.02500000]), ...  % Insertion of extensor digitorum longus (in foot: toe)
    norm([0.0032; -0.1381; 0.0276] - ... % Origin of extensor digitorum longus (in shank)
    [0.0289;-0.4007; 0.0072]), ... % Via point 1 of extensor digitorum longus (in shank)
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6
 
% Extensor hallucis longus
Model.Generic.Lcst(41,:) = [0, ... % No lenght embedded in segment 1
    norm([0.09700000; 0.03890000; -0.02110000] - ... % Via point 2 of extensor hallucis longus (in foot: calcanenous)
    [0.12930000; 0.03090000; -0.02570000]) + ... % Via point 3 of extensor hallucis longus (in foot: calcanenous)
    norm([0.12930000; 0.03090000; -0.02570000] - ... % Via point 3 of extensor hallucis longus (in foot: calcanenous)
    [0.17340000; 0.01390000; -0.02800000]) + ... % Via point 4 of extensor hallucis longus (in foot: calcanenous)
    norm([0.17340000; 0.01390000; -0.02800000] - ... % Via point 4 of extensor hallucis longus (in foot: calcanenous)
    (t_calcaneous2toe + [0.02980000; 0.00410000; -0.02450000])) + ... % Via point 5 of extensor hallucis longus (in foot: toe)
    norm([0.02980000; 0.00410000; -0.02450000] - ... % Via point 5 of extensor hallucis longus (in foot: toe)
    [0.05630000; 0.00340000; -0.01860000]), ... % Insertion of extensor hallucis longus (in foot: toe)
    norm([0.0012; -0.1767; 0.0228] - ... % Origin of extensor hallucis longus (in shank)
    [0.0326;-0.3985;-0.0085]),... % Via point 1 of extensor hallucis longus (in shank)
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Flexor digitorum longus
Model.Generic.Lcst(42,:) = [0, ... % No lenght embedded in segment 1
    norm([0.04360000; 0.03150000; -0.02800000] - ... % Via point 2 of flexor digitorum longus (in foot: calcaneous)
    [0.07080000; 0.01760000; -0.02630000]) + ... % Via point 3 of flexor digitorum longus (in foot: calcaneous)
    norm([0.07080000; 0.01760000; -0.02630000] - ... % Via point 3 of flexor digitorum longus (in foot: calcaneous)
    [0.16580000; -0.00810000; 0.01160000]) + ... % Via point 4 of flexor digitorum longus (in foot: calcaneous)
    norm([0.16580000; -0.00810000; 0.01160000] - ... % Via point 4 of flexor digitorum longus (in foot: calcaneous)
    (t_calcaneous2toe + [-0.00190000; -0.00780000; 0.01470000])) + ... % Via point 5 of flexor digitorum longus (in foot: toe)
    norm([-0.00190000; -0.00780000; 0.01470000] - ... % Via point 5 of flexor digitorum longus (in foot: toe)
    [0.02850000; -0.00710000; 0.02150000]) + ... % Via point 6 of flexor digitorum longus (in foot: toe)
    norm([0.02850000; -0.00710000; 0.02150000] - ... % Via point 6 of flexor digitorum longus (in foot: toe)
    [0.04410000; -0.00600000; 0.02420000]), ... % Insertion of flexor digitorum longus (in foot: toe)
    norm([-0.0083; -0.2046; -0.0018] - ... % Origin of flexor digitorum longus
    [-0.0154;-0.4051;-0.0196]), ... % Via point 1 of flexor digitorum longus
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6

% Flexor hallucis longus
Model.Generic.Lcst(43,:) = [0, ... % No lenght embedded in segment 1
    norm([0.03740000; 0.02760000; -0.02410000] - ... % Via point 2 of flexor hallucis longus (in foot: calcaneous)
    [0.10380000; 0.00680000; -0.02560000]) + ... % Via point 3 of flexor hallucis longus (in foot: calcaneous)
    norm([0.10380000; 0.00680000; -0.02560000] - ... % Via point 3 of flexor hallucis longus (in foot: calcaneous)
    [0.17260000; -0.00530000; -0.02690000]) + ... % Via point 4 of flexor hallucis longus (in foot: calcaneous)
    norm([0.17260000; -0.00530000; -0.02690000] - ... % Via point 4 of flexor hallucis longus (in foot: calcaneous)
    (t_calcaneous2toe + [0.01550000; -0.00640000; -0.02650000])) + ... % Via point 5 of flexor hallucis longus (in foot: toe)
    norm([0.01550000; -0.00640000; -0.02650000] - ... % Via point 5 of flexor hallucis longus (in foot: toe)
    [0.05620000; -0.01020000; -0.01810000]), ... % Origin of flexor hallucis longus (in foot: toe)
    norm([-0.0079; -0.2334; 0.0244] - ... % Origin of flexor hallucis longus
    [-0.0186;-0.4079;-0.0174]), ... % Via point of flexor hallucis longus
    0, 0, 0]; % No lenght embedded in segment 4, 5, or 6


% -------------------------------------------------------------------------
% Force max (Delp 1990 = Wickiewicz 1983 & Brand 1986)
% -------------------------------------------------------------------------
Model.Generic.Fmax(1,1) = 382;           % Gluteus maximus I
Model.Generic.Fmax(2,1) = 546;           % Gluteus maximus II
Model.Generic.Fmax(3,1) = 368;           % Gluteus maximus III
Model.Generic.Fmax(4,1) = 546;           % Gluteus medius I
Model.Generic.Fmax(5,1) = 383;           % Gluteus medius II
Model.Generic.Fmax(6,1) = 435;           % Gluteus medius III
Model.Generic.Fmax(7,1) = 180;           % Gluteus minimus I
Model.Generic.Fmax(8,1) = 190;           % Gluteus minimus II
Model.Generic.Fmax(9,1) = 215;           % Gluteus minimus III
Model.Generic.Fmax(10,1) = 418;          % Adductor longus
Model.Generic.Fmax(11,1) = 286;          % Adductor brevis
Model.Generic.Fmax(12,1) = 346;          % Adductor magnus I
Model.Generic.Fmax(13,1) = 312;          % Adductor magnus II
Model.Generic.Fmax(14,1) = 444;          % Adductor magnus III
Model.Generic.Fmax(15,1) = 177;          % Pectineus
Model.Generic.Fmax(16,1) = 429;          % Illiacus
Model.Generic.Fmax(17,1) = 371;          % Psoas
Model.Generic.Fmax(18,1) = 254;          % Quadratus femoris
Model.Generic.Fmax(19,1) = 109;          % Gemelli
Model.Generic.Fmax(20,1) = 296;          % Piriformis
Model.Generic.Fmax(21,1) = 155;          % Tensor fasciae latae
Model.Generic.Fmax(22,1) = 108;          % Gracilis
Model.Generic.Fmax(23,1) = 104;          % Sartorius
Model.Generic.Fmax(24,1) = 1030;         % Semimembranosus
Model.Generic.Fmax(25,1) = 328;          % Semitendinosus
Model.Generic.Fmax(26,1) = 717;          % Biceps femoris long head
Model.Generic.Fmax(27,1) = 402;          % Biceps femoris short head
Model.Generic.Fmax(28,1) = 779;          % Rectus femoris
Model.Generic.Fmax(29,1) = 1294;         % Vastus medialis
Model.Generic.Fmax(30,1) = 1365;         % Vastus intermedius
Model.Generic.Fmax(31,1) = 1871;         % Vastus lateralis
Model.Generic.Fmax(32,1) = 1113;         % Gastrocnemius medialis
Model.Generic.Fmax(33,1) = 488;          % Gastrocnemius lateralis
Model.Generic.Fmax(34,1) = 2839;         % Soleus
Model.Generic.Fmax(35,1) = 1270;         % Tibialis posterior
Model.Generic.Fmax(36,1) = 603;          % Tibialis anterior
Model.Generic.Fmax(37,1) = 348;          % Peroneus brevis
Model.Generic.Fmax(38,1) = 754;          % Peroneus longus
Model.Generic.Fmax(39,1) = 90;           % Peroneus tertius
Model.Generic.Fmax(40,1) = 341;          % Extensor digitorum longus
Model.Generic.Fmax(41,1) = 108;          % Extensor hallucis longus
Model.Generic.Fmax(42,1) = 310;          % Flexor digitorum longus
Model.Generic.Fmax(43,1) = 322;          % Flexor hallucis longus

% -------------------------------------------------------------------------
% PCSA (Thorpe 1997 + Friedrich and Brand 1990)
% -------------------------------------------------------------------------
Model.Generic.PCSA(1,1) = 41.7/3;        % Gluteus maximus I
Model.Generic.PCSA(2,1) = 41.7/3;        % Gluteus maximus II
Model.Generic.PCSA(3,1) = 41.7/3;        % Gluteus maximus III
Model.Generic.PCSA(4,1) = 21.60;         % Gluteus medius I
Model.Generic.PCSA(5,1) = 15.70;         % Gluteus medius II
Model.Generic.PCSA(6,1) = 16.80;         % Gluteus medius III
Model.Generic.PCSA(7,1) = 9.60;          % Gluteus minimus I
Model.Generic.PCSA(8,1) = 10.00;         % Gluteus minimus II
Model.Generic.PCSA(9,1) = 11.00;         % Gluteus minimus III
Model.Generic.PCSA(10,1) = 15.00;        % Adductor longus
Model.Generic.PCSA(11,1) = 9.70;         % Adductor brevis
Model.Generic.PCSA(12,1) = (22+25)/3;    % Adductor magnus I
Model.Generic.PCSA(13,1) = (22+25)/3;    % Adductor magnus II
Model.Generic.PCSA(14,1) = (22+25)/3;    % Adductor magnus III
Model.Generic.PCSA(15,1) = 5.40;         % Pectineus
Model.Generic.PCSA(16,1) = 32.00;        % Illiacus
Model.Generic.PCSA(17,1) = 15.6;         % Psoas
Model.Generic.PCSA(18,1) = 21.00;        % Quadratus femoris
Model.Generic.PCSA(19,1) = 4.33;         % Gemelli
Model.Generic.PCSA(20,1) = 20.54;        % Piriformis
Model.Generic.PCSA(21,1) = 5.20;         % Tensor fasciae latae
Model.Generic.PCSA(22,1) = 3.40;         % Gracilis
Model.Generic.PCSA(23,1) = 3.65;         % Sartorius
Model.Generic.PCSA(24,1) = 39.90;        % Semimembranosus
Model.Generic.PCSA(25,1) = 9.40;         % Semitendinosus
Model.Generic.PCSA(26,1) = 28.80;        % Biceps femoris long head
Model.Generic.PCSA(27,1) = 8.14;         % Biceps femoris short head
Model.Generic.PCSA(28,1) = 33.60;        % Rectus femoris
Model.Generic.PCSA(29,1) = 46.70;        % Vastus medialis
Model.Generic.PCSA(30,1) = 53.70;        % Vastus intermedius
Model.Generic.PCSA(31,1) = 68.80;        % Vastus lateralis
Model.Generic.PCSA(32,1) = 41.80;        % Gastrocnemius medialis
Model.Generic.PCSA(33,1) = 19.90;        % Gastrocnemius lateralis
Model.Generic.PCSA(34,1) = 118.70;       % Soleus
Model.Generic.PCSA(35,1) = 36.20;        % Tibialis posterior
Model.Generic.PCSA(36,1) = 20.40;        % Tibialis anterior
Model.Generic.PCSA(37,1) = 11.50;        % Peroneus brevis
Model.Generic.PCSA(38,1) = 21.40;        % Peroneus longus
Model.Generic.PCSA(39,1) = 4.14;         % Peroneus tertius
Model.Generic.PCSA(40,1) = 10.50;        % Extensor digitorum longus
Model.Generic.PCSA(41,1) = 4.85;         % Extensor hallucis longus
Model.Generic.PCSA(42,1) = 9.90;         % Flexor digitorum longus
Model.Generic.PCSA(43,1) = 14.10;        % Flexor hallucis longus

% -------------------------------------------------------------------------
% Optimal isometric length (m) (Delp 1990)
% -------------------------------------------------------------------------
Model.Generic.L0(1,1) = 0.14200000;     % Gluteus maximus I
Model.Generic.L0(2,1) = 0.14700000;     % Gluteus maximus II
Model.Generic.L0(3,1) = 0.14400000;     % Gluteus maximus III
Model.Generic.L0(4,1) = 0.05350000;     % Gluteus medius I
Model.Generic.L0(5,1) = 0.08450000;     % Gluteus medius II
Model.Generic.L0(6,1) = 0.06460000;     % Gluteus medius III
Model.Generic.L0(7,1) = 0.06800000;     % Gluteus minimus I
Model.Generic.L0(8,1) = 0.05600000;     % Gluteus minimus II
Model.Generic.L0(9,1) = 0.03800000;     % Gluteus minimus III
Model.Generic.L0(10,1) = 0.13800000;    % Adductor longus
Model.Generic.L0(11,1) = 0.13300000;    % Adductor brevis
Model.Generic.L0(12,1) = 0.08700000;    % Adductor magnus I
Model.Generic.L0(13,1) = 0.12100000;    % Adductor magnus II
Model.Generic.L0(14,1) = 0.13100000;    % Adductor magnus III
Model.Generic.L0(15,1) = 0.13300000;    % Pectineus
Model.Generic.L0(16,1) = 0.10000000;    % Iliacus
Model.Generic.L0(17,1) = 0.10400000;    % Psoas
Model.Generic.L0(18,1) = 0.05400000;    % Quadratus femoris
Model.Generic.L0(19,1) = 0.02400000;    % Gemelli
Model.Generic.L0(20,1) = 0.02600000;    % Piriformis
Model.Generic.L0(21,1) = 0.09500000;    % Tensor fascia latae
Model.Generic.L0(22,1) = 0.35200000;    % Gracilis
Model.Generic.L0(23,1) = 0.57900000;    % Sartorius
Model.Generic.L0(24,1) = 0.08000000;    % Semimembranosus
Model.Generic.L0(25,1) = 0.20100000;    % Semitendinosus
Model.Generic.L0(26,1) = 0.10900000;    % Biceps femoris long head
Model.Generic.L0(27,1) = 0.17300000;    % Biceps femoris short head
Model.Generic.L0(28,1) = 0.08400000;    % Rectus femoris
Model.Generic.L0(29,1) = 0.08900000;    % Vastus medialis
Model.Generic.L0(30,1) = 0.08700000;    % Vastus intermedialis
Model.Generic.L0(31,1) = 0.08400000;    % Vastus lateralis
Model.Generic.L0(32,1) = 0.04500000;    % Gastrocnemius medialis
Model.Generic.L0(33,1) = 0.06400000;    % Gastrocnemius lateralis
Model.Generic.L0(34,1) = 0.03000000;    % Soleus
Model.Generic.L0(35,1) = 0.03100000;    % Tibialis posterior
Model.Generic.L0(36,1) = 0.09800000;    % Tibialis anterior
Model.Generic.L0(37,1) = 0.05000000;    % Peroneus brevis
Model.Generic.L0(38,1) = 0.04900000;    % Peroneus longus
Model.Generic.L0(39,1) = 0.07900000;    % Peroneus tertius
Model.Generic.L0(40,1) = 0.10200000;    % Extensor digitorum longus
Model.Generic.L0(41,1) = 0.11100000;    % Extensor hallucis longus
Model.Generic.L0(42,1) = 0.03400000;    % Flexor digitorum longus
Model.Generic.L0(43,1) = 0.04300000;    % Flexor hallucis longus

% -------------------------------------------------------------------------
% Maximum contraction velocity of the fibers, in optimal fiberlengths
% per second (m/sec) (Delp 1990)
% -------------------------------------------------------------------------
for i = 1:43
    Model.Generic.V0max(i,1) = 10 * Model.Generic.L0(i,1);
end

% -------------------------------------------------------------------------
% Tendon slack length (m) (Delp 1990)
% -------------------------------------------------------------------------
Model.Generic.Lts(1,1) = 0.12500000;     % Gluteus maximus I
Model.Generic.Lts(2,1) = 0.12700000;     % Gluteus maximus II
Model.Generic.Lts(3,1) = 0.14500000;     % Gluteus maximus III
Model.Generic.Lts(4,1) = 0.07800000;     % Gluteus medius I
Model.Generic.Lts(5,1) = 0.05300000;     % Gluteus medius II
Model.Generic.Lts(6,1) = 0.05300000;     % Gluteus medius III
Model.Generic.Lts(7,1) = 0.01600000;     % Gluteus minimus I
Model.Generic.Lts(8,1) = 0.02600000;     % Gluteus minimus II
Model.Generic.Lts(9,1) = 0.05100000;     % Gluteus minimus III
Model.Generic.Lts(10,1) = 0.11000000;    % Adductor longus
Model.Generic.Lts(11,1) = 0.02000000;    % Adductor brevis
Model.Generic.Lts(12,1) = 0.06000000;    % Adductor magnus I
Model.Generic.Lts(13,1) = 0.13000000;    % Adductor magnus II
Model.Generic.Lts(14,1) = 0.26000000;    % Adductor magnus III
Model.Generic.Lts(15,1) = 0.00100000;    % Pectineus
Model.Generic.Lts(16,1) = 0.09000000;    % Iliacus
Model.Generic.Lts(17,1) = 0.13000000;    % Psoas
Model.Generic.Lts(18,1) = 0.02400000;    % Quadratus femoris
Model.Generic.Lts(19,1) = 0.03900000;    % Gemelli
Model.Generic.Lts(20,1) = 0.11500000;    % Piriformis
Model.Generic.Lts(21,1) = 0.42500000;    % Tensor fascia latae
Model.Generic.Lts(22,1) = 0.14000000;    % Gracilis
Model.Generic.Lts(23,1) = 0.04000000;    % Sartorius
Model.Generic.Lts(24,1) = 0.35900000;    % Semimembranosus
Model.Generic.Lts(25,1) = 0.26200000;    % Semitendinosus
Model.Generic.Lts(26,1) = 0.34100000;    % Biceps femoris long head
Model.Generic.Lts(27,1) = 0.10000000;    % Biceps femoris short head
Model.Generic.Lts(28,1) = 0.34600000;    % Rectus femoris
Model.Generic.Lts(29,1) = 0.12600000;    % Vastus medialis
Model.Generic.Lts(30,1) = 0.13600000;    % Vastus intermedialis
Model.Generic.Lts(31,1) = 0.15700000;    % Vastus lateralis
Model.Generic.Lts(32,1) = 0.40800000;    % Gastrocnemius medialis
Model.Generic.Lts(33,1) = 0.38500000;    % Gastrocnemius lateralis
Model.Generic.Lts(34,1) = 0.26800000;    % Soleus
Model.Generic.Lts(35,1) = 0.31000000;    % Tibialis posterior
Model.Generic.Lts(36,1) = 0.22300000;    % Tibialis anterior
Model.Generic.Lts(37,1) = 0.16100000;    % Peroneus brevis
Model.Generic.Lts(38,1) = 0.34500000;    % Peroneus longus
Model.Generic.Lts(39,1) = 0.10000000;    % Peroneus tertius
Model.Generic.Lts(40,1) = 0.34500000;    % Extensor digitorum longus
Model.Generic.Lts(41,1) = 0.30500000;    % Extensor hallucis longus
Model.Generic.Lts(42,1) = 0.40000000;    % Flexor digitorum longus
Model.Generic.Lts(43,1) = 0.38000000;    % Flexor hallucis longus

% -------------------------------------------------------------------------
% Pennation angle at optimal isometric length (radian) (Delp 1990)
% -------------------------------------------------------------------------
Model.Generic.pennation(1,1) = 0.08726646;           % Gluteus maximus I
Model.Generic.pennation(2,1) = 0.00000000;           % Gluteus maximus II
Model.Generic.pennation(3,1) = 0.08726646;           % Gluteus maximus III
Model.Generic.pennation(4,1) = 0.13962634;           % Gluteus medius I
Model.Generic.pennation(5,1) = 0.00000000;           % Gluteus medius II
Model.Generic.pennation(6,1) = 0.33161256;           % Gluteus medius III
Model.Generic.pennation(7,1) = 0.17453293;           % Gluteus minimus I
Model.Generic.pennation(8,1) = 0.00000000;           % Gluteus minimus II
Model.Generic.pennation(9,1) = 0.36651914;           % Gluteus minimus III
Model.Generic.pennation(10,1) = 0.10471976;          % Adductor longus
Model.Generic.pennation(11,1) = 0.00000000;          % Adductor brevis
Model.Generic.pennation(12,1) = 0.08726646;          % Adductor magnus I
Model.Generic.pennation(13,1) = 0.05235988;          % Adductor magnus II
Model.Generic.pennation(14,1) = 0.08726646;          % Adductor magnus III
Model.Generic.pennation(15,1) = 0.00000000;          % Pectineus
Model.Generic.pennation(16,1) = 0.12217305;          % Iliacus
Model.Generic.pennation(17,1) = 0.13962634;          % Psoas
Model.Generic.pennation(18,1) = 0.00000000;          % Quadratus femoris
Model.Generic.pennation(19,1) = 0.00000000;          % Gemelli
Model.Generic.pennation(20,1) = 0.17453293;          % Piriformis
Model.Generic.pennation(21,1) = 0.05235988;          % Tensor fascia latae
Model.Generic.pennation(22,1) = 0.05235988;          % Gracilis
Model.Generic.pennation(23,1) = 0.00000000;          % Sartorius
Model.Generic.pennation(24,1) = 0.26179939;          % Semimembranosus
Model.Generic.pennation(25,1) = 0.08726646;          % Semitendinosus
Model.Generic.pennation(26,1) = 0.00000000;          % Biceps femoris long head
Model.Generic.pennation(27,1) = 0.40142573;          % Biceps femoris short head
Model.Generic.pennation(28,1) = 0.08726646;          % Rectus femoris
Model.Generic.pennation(29,1) = 0.08726646;          % Vastus medialis
Model.Generic.pennation(30,1) = 0.05235988;          % Vastus intermedialis
Model.Generic.pennation(31,1) = 0.08726646;          % Vastus lateralis
Model.Generic.pennation(32,1) = 0.29670597;          % Gastrocnemius medialis
Model.Generic.pennation(33,1) = 0.13962634;          % Gastrocnemius lateralis
Model.Generic.pennation(34,1) = 0.43633231;          % Soleus
Model.Generic.pennation(35,1) = 0.20943951;          % Tibialis posterior
Model.Generic.pennation(36,1) = 0.08726646;          % Tibialis anterior
Model.Generic.pennation(37,1) = 0.08726646;          % Peroneus brevis
Model.Generic.pennation(38,1) = 0.17453293;          % Peroneus longus
Model.Generic.pennation(39,1) = 0.22689280;          % Peroneus tertius
Model.Generic.pennation(40,1) = 0.13962634;          % Extensor digitorum longus
Model.Generic.pennation(41,1) = 0.10471976;          % Extensor hallucis longus
Model.Generic.pennation(42,1) = 0.12217305;          % Flexor digitorum longus
Model.Generic.pennation(43,1) = 0.17453293;          % Flexor hallucis longus


%% -------------------------------------------------------------------------
% Bone models
% asc files from Delp
% -------------------------------------------------------------------------

% Foot bone
bone_file = "foot.asc";
fid = fopen(bone_file,'r');
fgetl(fid); % Header NORM_ASCII
nPoints = fscanf(fid,'%f',1);
nPolys = fscanf(fid,'%f',1);
fgetl(fid); % Go to next line
fgetl(fid); % Bounding box of the bone
for i = 1:nPoints
    Segment(2).Generic.Vertices(i,1) = fscanf(fid,'%f',1) + t_midmalleoli2calcaneus(1,1);
    Segment(2).Generic.Vertices(i,2) = fscanf(fid,'%f',1) + t_midmalleoli2calcaneus(2,1);
    Segment(2).Generic.Vertices(i,3) = fscanf(fid,'%f',1) + t_midmalleoli2calcaneus(3,1);
    fscanf(fid,'%f',1); 
    fscanf(fid,'%f',1);
    fscanf(fid,'%f',1); % Coordinates of the vertex normal
end
for i = 1:nPolys 
    nVect = fscanf(fid,'%f',1); % Face with 3 or 4 vertices
    Segment(2).Generic.Faces(i,1) = fscanf(fid,'%f',1)+1; % +1 to start count at 1 instead of 0
    Segment(2).Generic.Faces(i,2) = fscanf(fid,'%f',1)+1;
    Segment(2).Generic.Faces(i,3) = fscanf(fid,'%f',1)+1;
    if nVect == 3
        Segment(2).Generic.Faces(i,4) = NaN;
    else % nVect == 4
        Segment(2).Generic.Faces(i,4) = fscanf(fid,'%f',1)+1;
    end
end
fclose(fid);



% Tibial Insert
% From CT-scan/CAD femur origin to CT-scan/CAD tibia origin
Tfimplant2timplant = [eye(3),-[1.14; 41.09; 0]/1000; [0,0,0,1]];
%
[Triangles,Points] = stlread('Tibia_Fibula_And_Tibial_Insert.stl');
Points(:,3) = - Points(:,3); % Mesh definition (left knee)
Segment(3).Generic.Vertices = (inv(Segment(3).Generic.T_implant2bone)*Tfimplant2timplant*...
    [Points'/1000;ones(1,length(Points))])';
Segment(3).Generic.Vertices(:,4) = [];
Segment(3).Generic.Faces = Triangles;
            
% Patella bone
bone_file = "patella.asc";
fid = fopen(bone_file,'r');
fgetl(fid); % Header NORM_ASCII
nPoints = fscanf(fid,'%f',1);
nPolys = fscanf(fid,'%f',1);
fgetl(fid); % Go to next line
fgetl(fid); % Bounding box of the bone
for i = 1:nPoints
    Segment(4).Generic.Vertices(i,1) = fscanf(fid,'%f',1) + t_proximal2patella(1,1); 
    Segment(4).Generic.Vertices(i,2) = fscanf(fid,'%f',1) + t_proximal2patella(2,1);
    Segment(4).Generic.Vertices(i,3) = fscanf(fid,'%f',1) + t_proximal2patella(3,1);
    fscanf(fid,'%f',1); 
    fscanf(fid,'%f',1);
    fscanf(fid,'%f',1); % Coordinates of the vertex normal
end
for i = 1:nPolys 
    nVect = fscanf(fid,'%f',1); % Face with 3 or 4 vertices
    Segment(4).Generic.Faces(i,1) = fscanf(fid,'%f',1)+1; % +1 to start count at 1 instead of 0
    Segment(4).Generic.Faces(i,2) = fscanf(fid,'%f',1)+1;
    Segment(4).Generic.Faces(i,3) = fscanf(fid,'%f',1)+1;
    if nVect == 3
        Segment(4).Generic.Faces(i,4) = NaN;
    else % nVect == 4
        Segment(4).Generic.Faces(i,4) = fscanf(fid,'%f',1)+1;
    end
end
fclose(fid);

% Femoral component
[Triangles,Points] = stlread('Femur_And_Femoral_Component.stl');
Points(:,3) = - Points(:,3); % Mesh definition (left knee)
Segment(5).Generic.Vertices = (inv(Segment(5).Generic.T_implant2bone)*...
    [Points'/1000;ones(1,length(Points))])';
Segment(5).Generic.Vertices(:,4) = [];
Segment(5).Generic.Faces = Triangles;

% Pelvis bone
bone_file = "pelvis.asc";
fid = fopen(bone_file,'r');
fgetl(fid); % Header NORM_ASCII
nPoints = fscanf(fid,'%f',1);
nPolys = fscanf(fid,'%f',1);
fgetl(fid); % Go to next line
fgetl(fid); % Bounding box of the bone
for i = 1:nPoints
    Segment(6).Generic.Vertices(i,1) = fscanf(fid,'%f',1) + t_lumbar2pelvis(1,1);
    Segment(6).Generic.Vertices(i,2) = fscanf(fid,'%f',1) + t_lumbar2pelvis(2,1);
    Segment(6).Generic.Vertices(i,3) = fscanf(fid,'%f',1) + t_lumbar2pelvis(3,1);
    fscanf(fid,'%f',1); 
    fscanf(fid,'%f',1);
    fscanf(fid,'%f',1); % Coordinates of the vertex normal
end
for i = 1:nPolys 
    nVect = fscanf(fid,'%f',1); % Face with 3 or 4 vertices
    Segment(6).Generic.Faces(i,1) = fscanf(fid,'%f',1)+1; % +1 to start count at 1 instead of 0
    Segment(6).Generic.Faces(i,2) = fscanf(fid,'%f',1)+1;
    Segment(6).Generic.Faces(i,3) = fscanf(fid,'%f',1)+1;
    if nVect == 3
        Segment(6).Generic.Faces(i,4) = NaN;
    else % nVect == 4
        Segment(6).Generic.Faces(i,4) = fscanf(fid,'%f',1)+1;
    end
end
fclose(fid);


%% ------------------------------------------------------------------------
% Subject-specific segment length
% -------------------------------------------------------------------------

Segment(3).Generic.L = 0.4004; % From SC CT-scan
Segment(3).Generic.W = 0.0666; % From SC CT-scan
Segment(4).Generic.L = norm([45.22; -24.53; 0.05]/1000 - ... % Patellar tendon origin
    [39.0; -0.36; -2.79]/1000); % Origin of patella;
Segment(5).Generic.L = 0.4136; % From SC CT-scan
Segment(5).Generic.W = 0.0835; % From SC CT-scan
%
Joint(4).Generic.d(1,1) = 44.447/1000;


%% ------------------------------------------------------------------------
% Scale between Delp and DM
% -------------------------------------------------------------------------

Segment(3).Informed.Scale = Segment(3).Generic.L/...
    norm((LM(1:2,1) + MM(1:2,1))/2 - ... % From mid-epicondyles to mid-malleolus
    (LFE(1:2,1) + MFE(1:2,1))/2); % Default
Segment(5).Informed.Scale = Segment(5).Generic.L/...
    (norm([(LFE(1:2,1) + MFE(1:2,1))/2; 0] - ... % From hip joint centre to mid-epicondyles 
    HJC)); % Default
% Segment(5).Informed.Scale = Segment(3).Informed.Scale;

Segment(2).Informed.Scale = (Segment(3).Informed.Scale + ...
    Segment(5).Informed.Scale)/2; % Default
Segment(4).Informed.Scale = (Segment(3).Informed.Scale + ...
    Segment(5).Informed.Scale)/2; % Default
Segment(6).Informed.Scale = (Segment(3).Informed.Scale + ...
    Segment(5).Informed.Scale)/2; % Default

Segment(2).Generic.L = Segment(2).Generic.L*Segment(2).Informed.Scale;
Segment(2).Generic.W = Segment(2).Generic.W*Segment(2).Informed.Scale;
Segment(6).Generic.L = Segment(6).Generic.L*Segment(6).Informed.Scale;
Segment(6).Generic.W = Segment(6).Generic.W*Segment(6).Informed.Scale;

for i = [2,3,5,6]
    Segment(i).Generic.rMs = Segment(i).Generic.rMs*...
        Segment(i).Informed.Scale;
end
Segment(2).Generic.rVs = Segment(2).Generic.rVs*Segment(2).Informed.Scale;
Segment(3).Generic.rVs(:,[1,10:28]) = ... % Others are subject-specific
    Segment(3).Generic.rVs(:,[1,10:28])*Segment(3).Informed.Scale;
Segment(4).Generic.rVs(:,2:5) = ... % Other is subject-specific
    Segment(4).Generic.rVs(:,2:5)*Segment(4).Informed.Scale;
Segment(5).Generic.rVs(:,7:43) = ... % Others are subject-specific
    Segment(5).Generic.rVs(:,7:43)*Segment(5).Informed.Scale;
Segment(6).Generic.rVs = Segment(6).Generic.rVs*Segment(6).Informed.Scale;

