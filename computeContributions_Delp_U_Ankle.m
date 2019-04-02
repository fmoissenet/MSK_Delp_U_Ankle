function Contribution = computeContributions_Delp_U_Ankle(Segment,Joint,Model,weight)

% Number of frames
n = size(Segment(2).rM,3);

% Swing
temp = find(Joint(1).F(2,:,:)==0,2);
sw = temp(2);
clear temp;

% -------------------------------------------------------------------------
% Introduction of 3D forces
% -------------------------------------------------------------------------

% Prepare variables
% -------------------------------------------------------------------------
% Segment parameters
u1 = Segment(1).Q(1:3,1,:);
rP1 = Segment(1).Q(4:6,1,:);
rD1 = Segment(1).Q(7:9,1,:);
w1 = Segment(1).Q(10:12,1,:);

u2 = Segment(2).Q(1:3,1,:);
rP2 = Segment(2).Q(4:6,1,:);
rD2 = Segment(2).Q(7:9,1,:);
w2 = Segment(2).Q(10:12,1,:);

u3 = Segment(3).Q(1:3,1,:);
rP3 = Segment(3).Q(4:6,1,:);
rD3 = Segment(3).Q(7:9,1,:);
w3 = Segment(3).Q(10:12,1,:);

u4 = Segment(4).Q(1:3,1,:);
rP4 = Segment(4).Q(4:6,1,:);
rD4 = Segment(4).Q(7:9,1,:);
w4 = Segment(4).Q(10:12,1,:);

u5 = Segment(5).Q(1:3,1,:);
rP5 = Segment(5).Q(4:6,1,:);
rD5 = Segment(5).Q(7:9,1,:);
w5 = Segment(5).Q(10:12,1,:);

% Center of mass rCi
temp = Mprod_array3(Q2Tuw_array3(Segment(2).Q),...
    repmat([Segment(2).rCs;1],[1 1 n]));
rC2 = temp(1:3,:,:);
temp = Mprod_array3(Q2Tuw_array3(Segment(3).Q),...
    repmat([Segment(3).rCs;1],[1 1 n]));
rC3 = temp(1:3,:,:);
temp = Mprod_array3(Q2Tuw_array3(Segment(4).Q),...
    repmat([Segment(4).rCs;1],[1 1 n]));
rC4 = temp(1:3,:,:);
temp = Mprod_array3(Q2Tuw_array3(Segment(5).Q),...
    repmat([Segment(5).rCs;1],[1 1 n]));
rC5 = temp(1:3,:,:);

% Transpose of the interpolation matrix of rPi
NPit = [zeros(3,3,n);...
    repmat(eye(3,3),[1,1,n]);...
    zeros(3,3,n);...
    zeros(3,3,n)];

% Composition of interpolation matrices and lever arms
M33_2(1:3,2,:) = rP2-rD2;
M33_2(4:6,3,:) = -w2;
M33_2(7:9,3,:) = w2;
M33_2(10:12,1,:) = u2;

M33_3(1:3,2,:) = rP3-rD3;
M33_3(4:6,3,:) = -w3;
M33_3(7:9,3,:) = w3;
M33_3(10:12,1,:) = u3;

M33_4(1:3,2,:) = rP4-rD4;
M33_4(4:6,3,:) = -w4;
M33_4(7:9,3,:) = w4;
M33_4(10:12,1,:) = u4;

M33_5(1:3,2,:) = rP5-rD5;
M33_5(4:6,3,:) = -w5;
M33_5(7:9,3,:) = w5;
M33_5(10:12,1,:) = u5;

% Bstar
Bstar1(1:3,1,:) = cross(w1,u1);
Bstar1(1:3,2,:) = cross(u1,(rP1-rD1));
Bstar1(1:3,3,:) = cross(-(rP1-rD1),w1);

Bstar2(1:3,1,:) = cross(w2,u2);
Bstar2(1:3,2,:) = cross(u2,(rP2-rD2));
Bstar2(1:3,3,:) = cross(-(rP2-rD2),w2);

Bstar3(1:3,1,:) = cross(w3,u3);
Bstar3(1:3,2,:) = cross(u3,(rP3-rD3));
Bstar3(1:3,3,:) = cross(-(rP3-rD3),w3);

Bstar4(1:3,1,:) = cross(w4,u4);
Bstar4(1:3,2,:) = cross(u4,(rP4-rD4));
Bstar4(1:3,3,:) = cross(-(rP4-rD4),w4);

Bstar5(1:3,1,:) = cross(w5,u5);
Bstar5(1:3,2,:) = cross(u5,(rP5-rD5));
Bstar5(1:3,3,:) = cross(-(rP5-rD5),w5);

% Nstar2t (pseudo interpolation matrix)
Nstar2t = Mprod_array3(M33_2,Minv_array3(Bstar2));
Nstar3t = Mprod_array3(M33_3,Minv_array3(Bstar3));
Nstar4t = Mprod_array3(M33_4,Minv_array3(Bstar4));
Nstar5t = Mprod_array3(M33_5,Minv_array3(Bstar5));

% Transformation of ground reaction wrench {F0R 
% at center of pressure rP1                {M0R
% -------------------------------------------------------------------------
% Initialisation
Model.LR = zeros(48,6,n);

% LR
Model.LR(1:12,1:6,:) = [NPit + ...
    Mprod_array3(Nstar2t,Vskew_array3(rP1-rP2)), ...
    Mprod_array3(Nstar2t,Bstar1)];

% Transformation of dynamic wrenches {FiG 
% at center of mass rci              {MiG
% -------------------------------------------------------------------------
% Initialisation
Model.LG = zeros(48,24,n);

% LG
Model.LG(1:12,1:6,:) = [NPit + ...
    Mprod_array3(Nstar2t,Vskew_array3(rC2-rP2)), ...
    Mprod_array3(Nstar2t,Bstar2)];
Model.LG(13:24,7:12,:) = [NPit + ...
    Mprod_array3(Nstar3t,Vskew_array3(rC3-rP3)), ...
    Mprod_array3(Nstar3t,Bstar3)];
Model.LG(25:36,13:18,:) = [NPit + ...
    Mprod_array3(Nstar4t,Vskew_array3(rC4-rP4)), ...
    Mprod_array3(Nstar4t,Bstar4)];
Model.LG(37:48,19:24,:) = [NPit + ...
    Mprod_array3(Nstar5t,Vskew_array3(rC5-rP5)), ...
    Mprod_array3(Nstar5t,Bstar5)];

clear temp Contribution

% Compute contribution of weight to ground reaction and accelerations
% -------------------------------------------------------------------------
A = zeros(5,30,n);
b = zeros(5,1,n);
Astar = zeros(5,30,n);

A(:,:,1:sw-1) = [-Mprod_array3(Mtran_array3(Model.ZK(:,:,1:sw-1)),Model.LR(:,:,1:sw-1)) ...
    Mprod_array3(Mtran_array3(Model.ZK(:,:,1:sw-1)),Model.LG(:,:,1:sw-1))]; % Stance
A(:,:,sw:n) = [zeros(5,6,n-sw+1) ...
    Mprod_array3(Mtran_array3(Model.ZK(:,:,sw:n)),Model.LG(:,:,sw:n))]; % Swing
b(:,:,1:n) = Mprod_array3(Mtran_array3(Model.ZK),Model.P);

w = 1e-2; % Scale factor
Q0 = repmat(eye(30),[1,1,n]);
Q0(1,1,:) = repmat(w,[1,1,n]);
Q0(2,2,:) = repmat(w,[1,1,n]);
Q0(3,3,:) = repmat(w,[1,1,n]);
Q0(4,4,:) = repmat(w,[1,1,n]);
Q0(5,5,:) = repmat(w,[1,1,n]);

Astar = Mprod_array3(A,Minv_array3(Q0));
temp = Mprod_array3(Mprod_array3(Minv_array3(Q0),Mpinv_array3(Astar)),b);

Contribution.weight.F1R = temp(1:3,:,:);
Contribution.weight.M1R = Mprod_array3(Bstar1,temp(4:6,:,:));    
Contribution.weight.d2Qdt2 = Mprod_array3(...
    Mpinv_array3(Model.G),Mprod_array3(Model.LG,temp(7:30,:,:)));

% Hip adductors
i = [1:25,65:66];
Contribution.muscleForce.O_hipadductors = [];
Contribution.muscleForce.O_hipadductors = contributionAnalysis(Segment,Contribution.muscleForce.O_hipadductors,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);

% Hip abductors
i = [50:61,62:64,81,102:103];
Contribution.muscleForce.O_hipabductors = [];
Contribution.muscleForce.O_hipabductors = contributionAnalysis(Segment,Contribution.muscleForce.O_hipabductors,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);

% Hip extensors
i = [38:49];
Contribution.muscleForce.O_hipextensors = [];
Contribution.muscleForce.O_hipextensors = contributionAnalysis(Segment,Contribution.muscleForce.O_hipextensors,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);

% Hip flexors
i = [67,84:85];
Contribution.muscleForce.O_hipflexors = [];
Contribution.muscleForce.O_hipflexors = contributionAnalysis(Segment,Contribution.muscleForce.O_hipflexors,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);

% Knee extensors
i = [90:91,106:129];
Contribution.muscleForce.O_kneeextensors = [];
Contribution.muscleForce.O_kneeextensors = contributionAnalysis(Segment,Contribution.muscleForce.O_kneeextensors,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);

% Knee flexors
i = [26:29,94:95];
Contribution.muscleForce.O_kneeflexors = [];
Contribution.muscleForce.O_kneeflexors = contributionAnalysis(Segment,Contribution.muscleForce.O_kneeflexors,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);

% Ankle plantarflexors
i = [34:35,96:101];
Contribution.muscleForce.O_ankleplantarflex = [];
Contribution.muscleForce.O_ankleplantarflex = contributionAnalysis(Segment,Contribution.muscleForce.O_ankleplantarflex,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);

% Ankle dorsiflexors
i = [31:32,104];
Contribution.muscleForce.O_ankledorsiflex = [];
Contribution.muscleForce.O_ankledorsiflex = contributionAnalysis(Segment,Contribution.muscleForce.O_ankledorsiflex,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);

% Ankle eversors
i = [78:80];
Contribution.muscleForce.O_ankleeversors = [];
Contribution.muscleForce.O_ankleeversors = contributionAnalysis(Segment,Contribution.muscleForce.O_ankleeversors,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);

% Ankle inversors
i = [105];
Contribution.muscleForce.O_ankleinversors = [];
Contribution.muscleForce.O_ankleinversors = contributionAnalysis(Segment,Contribution.muscleForce.O_ankleinversors,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductors FULL
% i = 1:25;
% Contribution.muscleForce.adductorsFULL = [];
% Contribution.muscleForce.adductorsFULL = contributionAnalysis(Segment,Contribution.muscleForce.adductorsFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductor brevis FULL
% i = 1:6;  
% Contribution.muscleForce.adductorbrevisFULL = [];
% Contribution.muscleForce.adductorbrevisFULL = contributionAnalysis(Segment,Contribution.muscleForce.adductorbrevisFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductor brevis prox
% i = 1:2;  
% Contribution.muscleForce.adductorbrevisprox = [];
% Contribution.muscleForce.adductorbrevisprox = contributionAnalysis(Segment,Contribution.muscleForce.adductorbrevisprox,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductor brevis mid
% i = 3:4;  
% Contribution.muscleForce.adductorbrevismid = [];
% Contribution.muscleForce.adductorbrevismid = contributionAnalysis(Segment,Contribution.muscleForce.adductorbrevismid,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductor brevis dist
% i = 5:6;  
% Contribution.muscleForce.adductorbrevisdist = [];
% Contribution.muscleForce.adductorbrevisdist = contributionAnalysis(Segment,Contribution.muscleForce.adductorbrevisdist,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductor longus
% i = 7:12;  
% Contribution.muscleForce.adductorlongus = [];
% Contribution.muscleForce.adductorlongus = contributionAnalysis(Segment,Contribution.muscleForce.adductorlongus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductor magnus FULL
% i = 13:25;  
% Contribution.muscleForce.adductormagnusFULL = [];
% Contribution.muscleForce.adductormagnusFULL = contributionAnalysis(Segment,Contribution.muscleForce.adductormagnusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductor magnus dist
% i = 13:15;  
% Contribution.muscleForce.adductormagnusdist = [];
% Contribution.muscleForce.adductormagnusdist = contributionAnalysis(Segment,Contribution.muscleForce.adductormagnusdist,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductor magnus mid
% i = 16:21;  
% Contribution.muscleForce.adductormagnusmid = [];
% Contribution.muscleForce.adductormagnusmid = contributionAnalysis(Segment,Contribution.muscleForce.adductormagnusmid,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Adductor magnus prox
% i = 22:25;  
% Contribution.muscleForce.adductormagnusprox = [];
% Contribution.muscleForce.adductormagnusprox = contributionAnalysis(Segment,Contribution.muscleForce.adductormagnusprox,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Biceps femoris CL
% i = 26;  
% Contribution.muscleForce.bicepsfemorisCL = [];
% Contribution.muscleForce.bicepsfemorisCL = contributionAnalysis(Segment,Contribution.muscleForce.bicepsfemorisCL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Biceps femoris CB
% i = 27:29;  
% Contribution.muscleForce.bicepsfemorisCB = [];
% Contribution.muscleForce.bicepsfemorisCB = contributionAnalysis(Segment,Contribution.muscleForce.bicepsfemorisCB,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Extensor digitorum longus
% i = 30;  
% Contribution.muscleForce.extensordigitorumlongus = [];
% Contribution.muscleForce.extensordigitorumlongus = contributionAnalysis(Segment,Contribution.muscleForce.extensordigitorumlongus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Extensor hallux longus
% i = 31;  
% Contribution.muscleForce.extensorhalluxlongus = [];
% Contribution.muscleForce.extensorhalluxlongus = contributionAnalysis(Segment,Contribution.muscleForce.extensorhalluxlongus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Flexor digitorum longus
% i = 32;  
% Contribution.muscleForce.flexordigitorumlongus = [];
% Contribution.muscleForce.flexordigitorumlongus = contributionAnalysis(Segment,Contribution.muscleForce.flexordigitorumlongus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Flexor hallux longus
% i = 33;  
% Contribution.muscleForce.flexorhalluxlongus = [];
% Contribution.muscleForce.flexorhalluxlongus = contributionAnalysis(Segment,Contribution.muscleForce.flexorhalluxlongus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gastrocnemius FULL
% i = 34:35;  
% Contribution.muscleForce.gastrocnemiusFULL = [];
% Contribution.muscleForce.gastrocnemiusFULL = contributionAnalysis(Segment,Contribution.muscleForce.gastrocnemiusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gastrocnemius lateralis
% i = 34;  
% Contribution.muscleForce.gastrocnemiuslateralis = [];
% Contribution.muscleForce.gastrocnemiuslateralis = contributionAnalysis(Segment,Contribution.muscleForce.gastrocnemiuslateralis,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gastrocnemius medialis
% i = 35;  
% Contribution.muscleForce.gastrocnemiusmedialis = [];
% Contribution.muscleForce.gastrocnemiusmedialis = contributionAnalysis(Segment,Contribution.muscleForce.gastrocnemiusmedialis,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gemellus FULL
% i = 36:37;  
% Contribution.muscleForce.gemellusFULL = [];
% Contribution.muscleForce.gemellusFULL = contributionAnalysis(Segment,Contribution.muscleForce.gemellusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gemellus inf
% i = 36;  
% Contribution.muscleForce.gemellusinf = [];
% Contribution.muscleForce.gemellusinf = contributionAnalysis(Segment,Contribution.muscleForce.gemellusinf,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gemellus sup
% i = 37;  
% Contribution.muscleForce.gemellussup = [];
% Contribution.muscleForce.gemellussup = contributionAnalysis(Segment,Contribution.muscleForce.gemellussup,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus maximus FULL
% i = 38:49;  
% Contribution.muscleForce.gluteusmaximusFULL = [];
% Contribution.muscleForce.gluteusmaximusFULL = contributionAnalysis(Segment,Contribution.muscleForce.gluteusmaximusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus maximus sup
% i = 38:43;  
% Contribution.muscleForce.gluteusmaximussup = [];
% Contribution.muscleForce.gluteusmaximussup = contributionAnalysis(Segment,Contribution.muscleForce.gluteusmaximussup,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus maximus inf
% i = 44:49;  
% Contribution.muscleForce.gluteusmaximusinf = [];
% Contribution.muscleForce.gluteusmaximusinf = contributionAnalysis(Segment,Contribution.muscleForce.gluteusmaximusinf,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus medius FULL
% i = 50:61;  
% Contribution.muscleForce.gluteusmediusFULL = [];
% Contribution.muscleForce.gluteusmediusFULL = contributionAnalysis(Segment,Contribution.muscleForce.gluteusmediusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus medius ant
% i = 50:55;  
% Contribution.muscleForce.gluteusmediusant = [];
% Contribution.muscleForce.gluteusmediusant = contributionAnalysis(Segment,Contribution.muscleForce.gluteusmediusant,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus medius post
% i = 56:61;  
% Contribution.muscleForce.gluteusmediuspost = [];
% Contribution.muscleForce.gluteusmediuspost = contributionAnalysis(Segment,Contribution.muscleForce.gluteusmediuspost,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus minimus FULL
% i = 62:64;  
% Contribution.muscleForce.gluteusminimusFULL = [];
% Contribution.muscleForce.gluteusminimusFULL = contributionAnalysis(Segment,Contribution.muscleForce.gluteusminimusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus minimus lat
% i = 62;  
% Contribution.muscleForce.gluteusminimuslat = [];
% Contribution.muscleForce.gluteusminimuslat = contributionAnalysis(Segment,Contribution.muscleForce.gluteusminimuslat,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus minimus mid
% i = 63;  
% Contribution.muscleForce.gluteusminimusmid = [];
% Contribution.muscleForce.gluteusminimusmid = contributionAnalysis(Segment,Contribution.muscleForce.gluteusminimusmid,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gluteus minimus med
% i = 64;  
% Contribution.muscleForce.gluteusminimusmed = [];
% Contribution.muscleForce.gluteusminimusmed = contributionAnalysis(Segment,Contribution.muscleForce.gluteusminimusmed,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Gracilis
% i = 65:66;  
% Contribution.muscleForce.gracilis = [];
% Contribution.muscleForce.gracilis = contributionAnalysis(Segment,Contribution.muscleForce.gracilis,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Iliacus
% i = 67;  
% Contribution.muscleForce.iliacus = [];
% Contribution.muscleForce.iliacus = contributionAnalysis(Segment,Contribution.muscleForce.iliacus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Obturator externus FULL
% i = 68:72;  
% Contribution.muscleForce.obturatorexternusFULL = [];
% Contribution.muscleForce.obturatorexternusFULL = contributionAnalysis(Segment,Contribution.muscleForce.obturatorexternusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Obturator externus inf
% i = 68:69;  
% Contribution.muscleForce.obturatorexternusinf = [];
% Contribution.muscleForce.obturatorexternusinf = contributionAnalysis(Segment,Contribution.muscleForce.obturatorexternusinf,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Obturator externus sup
% i = 70:72;  
% Contribution.muscleForce.obturatorexternussup = [];
% Contribution.muscleForce.obturatorexternussup = contributionAnalysis(Segment,Contribution.muscleForce.obturatorexternussup,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Obturator internus
% i = 73;  
% Contribution.muscleForce.obturatorinternus = [];
% Contribution.muscleForce.obturatorinternus = contributionAnalysis(Segment,Contribution.muscleForce.obturatorinternus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Pectinus
% i = 74:77;  
% Contribution.muscleForce.pectinus = [];
% Contribution.muscleForce.pectinus = contributionAnalysis(Segment,Contribution.muscleForce.pectinus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Peroneus FULL
% i = 78:80;  
% Contribution.muscleForce.peroneusFULL = [];
% Contribution.muscleForce.peroneusFULL = contributionAnalysis(Segment,Contribution.muscleForce.peroneusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Peroneus brevis
% i = 78;  
% Contribution.muscleForce.peroneusbrevis = [];
% Contribution.muscleForce.peroneusbrevis = contributionAnalysis(Segment,Contribution.muscleForce.peroneusbrevis,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Peroneus longus
% i = 79;  
% Contribution.muscleForce.peroneuslongus = [];
% Contribution.muscleForce.peroneuslongus = contributionAnalysis(Segment,Contribution.muscleForce.peroneuslongus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Peroneus tertius
% i = 80;  
% Contribution.muscleForce.peroneustertius = [];
% Contribution.muscleForce.peroneustertius = contributionAnalysis(Segment,Contribution.muscleForce.peroneustertius,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Piriformis
% i = 81;  
% Contribution.muscleForce.piriformis = [];
% Contribution.muscleForce.piriformis = contributionAnalysis(Segment,Contribution.muscleForce.piriformis,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Plantaris
% i = 82;  
% Contribution.muscleForce.plantaris = [];
% Contribution.muscleForce.plantaris = contributionAnalysis(Segment,Contribution.muscleForce.plantaris,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Popliteus
% i = 83;  
% Contribution.muscleForce.popliteus = [];
% Contribution.muscleForce.popliteus = contributionAnalysis(Segment,Contribution.muscleForce.popliteus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Psoas minor
% i = 84;  
% Contribution.muscleForce.psoasminor = [];
% Contribution.muscleForce.psoasminor = contributionAnalysis(Segment,Contribution.muscleForce.psoasminor,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Iliopsoas FULL (psoas major + iliacus)
% i = [67;85];  
% Contribution.muscleForce.iliopsoasFULL = [];
% Contribution.muscleForce.iliopsoasFULL = contributionAnalysis(Segment,Contribution.muscleForce.iliopsoasFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Psoas major
% i = 85;  
% Contribution.muscleForce.psoasmajor = [];
% Contribution.muscleForce.psoasmajor = contributionAnalysis(Segment,Contribution.muscleForce.psoasmajor,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Quadratis femoris
% i = 86:89;  
% Contribution.muscleForce.quadratisfemoris = [];
% Contribution.muscleForce.quadratisfemoris = contributionAnalysis(Segment,Contribution.muscleForce.quadratisfemoris,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Rectus femoris
% i = 90:91;  
% Contribution.muscleForce.rectusfemoris = [];
% Contribution.muscleForce.rectusfemoris = contributionAnalysis(Segment,Contribution.muscleForce.rectusfemoris,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Sartorius FULL
% i = 92:93;  
% Contribution.muscleForce.sartoriusFULL = [];
% Contribution.muscleForce.sartoriusFULL = contributionAnalysis(Segment,Contribution.muscleForce.sartoriusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Sartorius prox
% i = 92;  
% Contribution.muscleForce.sartoriusprox = [];
% Contribution.muscleForce.sartoriusprox = contributionAnalysis(Segment,Contribution.muscleForce.sartoriusprox,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Sartorius dist
% i = 93;  
% Contribution.muscleForce.sartoriusdist = [];
% Contribution.muscleForce.sartoriusdist = contributionAnalysis(Segment,Contribution.muscleForce.sartoriusdist,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Hamstrings FULL
% i = [26:29 94:95];  
% Contribution.muscleForce.hamstringsFULL = [];
% Contribution.muscleForce.hamstringsFULL = contributionAnalysis(Segment,Contribution.muscleForce.hamstringsFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Semimembranosus
% i = 94;  
% Contribution.muscleForce.semimembranosus = [];
% Contribution.muscleForce.semimembranosus = contributionAnalysis(Segment,Contribution.muscleForce.semimembranosus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Semitendinosus
% i = 95;  
% Contribution.muscleForce.semitendinosus = [];
% Contribution.muscleForce.semitendinosus = contributionAnalysis(Segment,Contribution.muscleForce.semitendinosus,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Soleus FULL
% i = 96:101;  
% Contribution.muscleForce.soleusFULL = [];
% Contribution.muscleForce.soleusFULL = contributionAnalysis(Segment,Contribution.muscleForce.soleusFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Soleus med
% i = 96:98;  
% Contribution.muscleForce.soleusmed = [];
% Contribution.muscleForce.soleusmed = contributionAnalysis(Segment,Contribution.muscleForce.soleusmed,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Soleus lat
% i = 99:101;  
% Contribution.muscleForce.soleuslat = [];
% Contribution.muscleForce.soleuslat = contributionAnalysis(Segment,Contribution.muscleForce.soleuslat,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Tensor fascia lata
% i = 102:103;  
% Contribution.muscleForce.tensorfascialata = [];
% Contribution.muscleForce.tensorfascialata = contributionAnalysis(Segment,Contribution.muscleForce.tensorfascialata,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Tibialis anterior
% i = 104;  
% Contribution.muscleForce.tibialisanterior = [];
% Contribution.muscleForce.tibialisanterior = contributionAnalysis(Segment,Contribution.muscleForce.tibialisanterior,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Tibialis posterior
% i = 105;  
% Contribution.muscleForce.tibialisposterior = [];
% Contribution.muscleForce.tibialisposterior = contributionAnalysis(Segment,Contribution.muscleForce.tibialisposterior,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Vastii FULL
% i = 106:129;  
% Contribution.muscleForce.vastiiFULL = [];
% Contribution.muscleForce.vastiiFULL = contributionAnalysis(Segment,Contribution.muscleForce.vastiiFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Vastus intermedius
% i = 106:111;  
% Contribution.muscleForce.vastusintermedius = [];
% Contribution.muscleForce.vastusintermedius = contributionAnalysis(Segment,Contribution.muscleForce.vastusintermedius,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Vastus lateralis FULL
% i = 112:119;  
% Contribution.muscleForce.vastuslateralisFULL = [];
% Contribution.muscleForce.vastuslateralisFULL = contributionAnalysis(Segment,Contribution.muscleForce.vastuslateralisFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Vastus lateralis inf
% i = 112:117;  
% Contribution.muscleForce.vastuslateralisinf = [];
% Contribution.muscleForce.vastuslateralisinf = contributionAnalysis(Segment,Contribution.muscleForce.vastuslateralisinf,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Vastus lateralis sup
% i = 118:119;  
% Contribution.muscleForce.vastuslateralissup = [];
% Contribution.muscleForce.vastuslateralissup = contributionAnalysis(Segment,Contribution.muscleForce.vastuslateralissup,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Vastus medialis FULL
% i = 120:129;  
% Contribution.muscleForce.vastusmedialisFULL = [];
% Contribution.muscleForce.vastusmedialisFULL = contributionAnalysis(Segment,Contribution.muscleForce.vastusmedialisFULL,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Vastus medialis inf
% i = 120:121;  
% Contribution.muscleForce.vastusmedialisinf = [];
% Contribution.muscleForce.vastusmedialisinf = contributionAnalysis(Segment,Contribution.muscleForce.vastusmedialisinf,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Vastus medialis mid
% i = 122:123;  
% Contribution.muscleForce.vastusmedialismid = [];
% Contribution.muscleForce.vastusmedialismid = contributionAnalysis(Segment,Contribution.muscleForce.vastusmedialismid,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);
% 
% % Vastus medialis sup
% i = 124:129;  
% Contribution.muscleForce.vastusmedialissup = [];
% Contribution.muscleForce.vastusmedialissup = contributionAnalysis(Segment,Contribution.muscleForce.vastusmedialissup,Model,Bstar1,NPit,Nstar2t,rP1,rP2,i,sw,n);