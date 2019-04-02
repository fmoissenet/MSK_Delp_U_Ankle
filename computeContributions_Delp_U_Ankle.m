function Contribution = computeContributions_Delp_U_Ankle(Segment,Joint,Model,weight)

% Number of frames
n = size(Segment(2).rM,3);

% Swing (force threshold = 5 N)
temp = find(abs(Joint(1).F(2,:,10:end))<5,2);
sw   = temp(2);
clear temp;

% -------------------------------------------------------------------------
% Introduction of 3D forces
% -------------------------------------------------------------------------

% Prepare variables
% -------------------------------------------------------------------------
% Segment parameters
u1  = Segment(1).Q(1:3,1,:);
rP1 = Segment(1).Q(4:6,1,:);
rD1 = Segment(1).Q(7:9,1,:);
w1  = Segment(1).Q(10:12,1,:);

u2  = Segment(2).Q(1:3,1,:);
rP2 = Segment(2).Q(4:6,1,:);
rD2 = Segment(2).Q(7:9,1,:);
w2  = Segment(2).Q(10:12,1,:);

u3  = Segment(3).Q(1:3,1,:);
rP3 = Segment(3).Q(4:6,1,:);
rD3 = Segment(3).Q(7:9,1,:);
w3  = Segment(3).Q(10:12,1,:);

u4  = Segment(4).Q(1:3,1,:);
rP4 = Segment(4).Q(4:6,1,:);
rD4 = Segment(4).Q(7:9,1,:);
w4  = Segment(4).Q(10:12,1,:);

u5  = Segment(5).Q(1:3,1,:);
rP5 = Segment(5).Q(4:6,1,:);
rD5 = Segment(5).Q(7:9,1,:);
w5  = Segment(5).Q(10:12,1,:);

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
clear temp;

% Transpose of the interpolation matrix of rPi
NPit = [zeros(3,3,n);...
    repmat(eye(3,3),[1,1,n]);...
    zeros(3,3,n);...
    zeros(3,3,n)];

% Composition of interpolation matrices and lever arms
M33_2(1:3,2,:)   = rP2-rD2;
M33_2(4:6,3,:)   = -w2;
M33_2(7:9,3,:)   = w2;
M33_2(10:12,1,:) = u2;

M33_3(1:3,2,:)   = rP3-rD3;
M33_3(4:6,3,:)   = -w3;
M33_3(7:9,3,:)   = w3;
M33_3(10:12,1,:) = u3;

M33_4(1:3,2,:)   = rP4-rD4;
M33_4(4:6,3,:)   = -w4;
M33_4(7:9,3,:)   = w4;
M33_4(10:12,1,:) = u4;

M33_5(1:3,2,:)   = rP5-rD5;
M33_5(4:6,3,:)   = -w5;
M33_5(7:9,3,:)   = w5;
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
Model.LG(1:12,1:6,:)    = [NPit + ...
                           Mprod_array3(Nstar2t,Vskew_array3(rC2-rP2)), ...
                           Mprod_array3(Nstar2t,Bstar2)];
Model.LG(13:24,7:12,:)  = [NPit + ...
                           Mprod_array3(Nstar3t,Vskew_array3(rC3-rP3)), ...
                           Mprod_array3(Nstar3t,Bstar3)];
Model.LG(25:36,13:18,:) = [NPit + ...
                           Mprod_array3(Nstar4t,Vskew_array3(rC4-rP4)), ...
                           Mprod_array3(Nstar4t,Bstar4)];
Model.LG(37:48,19:24,:) = [NPit + ...
                           Mprod_array3(Nstar5t,Vskew_array3(rC5-rP5)), ...
                           Mprod_array3(Nstar5t,Bstar5)];

% -------------------------------------------------------------------------
% Compute contributions
% -------------------------------------------------------------------------                       
list ={'weight', ...
       'O_hipadductors', ...
       'O_hipabductors', ...
       'O_hipextensors', ...
       'O_hipflexors', ...
       'O_kneeextensors', ...
       'O_kneeflexors', ...
       'O_ankleplantarflex', ...
       'O_ankledorsiflex', ...
       'O_ankleeversors', ...
       'O_ankleinversors'};                       
                       
% Compute contribution of weight to ground reaction and accelerations
% -------------------------------------------------------------------------
A             = zeros(5,30,n);
b             = zeros(5,1,n);
Astar         = zeros(5,30,n);
A(:,:,2:sw-1) = [-Mprod_array3(Mtran_array3(Model.ZK(:,:,2:sw-1)),Model.LR(:,:,2:sw-1)) ...
                  Mprod_array3(Mtran_array3(Model.ZK(:,:,2:sw-1)),Model.LG(:,:,2:sw-1))]; % Stance
A(:,:,sw:n)   = [ zeros(5,6,n-sw+1) ...
                  Mprod_array3(Mtran_array3(Model.ZK(:,:,sw:n)),Model.LG(:,:,sw:n))]; % Swing
b(:,:,1:n)    = Mprod_array3(Mtran_array3(Model.ZK),Model.P);

% Weighted pseudo-inverse
w         = 1e-2; % Scale factor
Q0        = repmat(eye(30),[1,1,n]);
Q0(1,1,:) = repmat(w,[1,1,n]);
Q0(2,2,:) = repmat(w,[1,1,n]);
Q0(3,3,:) = repmat(w,[1,1,n]);
Q0(4,4,:) = repmat(w,[1,1,n]);
Q0(5,5,:) = repmat(w,[1,1,n]);
Astar     = Mprod_array3(A,Minv_array3(Q0));
temp      = Mprod_array3(Mprod_array3(Minv_array3(Q0),Mpinv_array3(Astar)),b);

% Contribution to ground reaction forces and moments and accelerations
Contribution.weight.F1R    = temp(1:3,:,:);
Contribution.weight.M1R    = Mprod_array3(Bstar1,temp(4:6,:,:));    
Contribution.weight.d2Qdt2 = Mprod_array3(Mpinv_array3(Model.G),...
                                          Mprod_array3(Model.LG,...
                                                       temp(7:30,:,:)));

% Compute contribution of musculo-tendon forces to ground reaction and acc.
% -------------------------------------------------------------------------                                                                                           
for m = 1:10
    if m == 1
        i = [10:14,22]; % Hip adductors
    elseif m == 2
        i = [4:6,7:9,20,21]; % Hip abductors
    elseif m == 3
        i = [1:3]; % Hip extensors
    elseif m == 4
        i = [16,17]; % Hip flexors
    elseif m == 5
        i = [28:31]; % Knee extensors
    elseif m == 6
        i = [24:27]; % Knee flexors
    elseif m == 7
        i = [32:34]; % Ankle plantarflexors
    elseif m == 8
        i = [40,41,36]; % Ankle dorsiflexors
    elseif m == 9
        i = [37:39]; % Ankle eversors
    elseif m == 10
        i = [35]; % Ankle inversors
    end
    A             = zeros(5,30,n);
    b             = zeros(5,1,n);
    Astar         = zeros(5,30,n);    
    A(:,:,2:sw-1) = [-Mprod_array3(Mtran_array3(Model.ZK(:,:,2:sw-1)),Model.LR(:,:,2:sw-1)) ...
        Mprod_array3(Mtran_array3(Model.ZK(:,:,2:sw-1)),Model.LG(:,:,2:sw-1))]; % Stance
    A(:,:,sw:n)   = [zeros(5,6,n-sw+1) ...
        Mprod_array3(Mtran_array3(Model.ZK(:,:,sw:n)),Model.LG(:,:,sw:n))]; % Swing
    b(:,:,1:n)    = Mprod_array3(Mtran_array3(Model.ZK),Mprod_array3(Model.Lever(:,i,:),Model.X(i,:,:)));
    
    % Weighted pseudo-inverse
    w         = 1e-2; % Scale factor
    Q0        = repmat(eye(30),[1,1,n]);
    Q0(1,1,:) = repmat(w,[1,1,n]);
    Q0(2,2,:) = repmat(w,[1,1,n]);
    Q0(3,3,:) = repmat(w,[1,1,n]);
    Q0(4,4,:) = repmat(w,[1,1,n]);
    Q0(5,5,:) = repmat(w,[1,1,n]);
    Astar     = Mprod_array3(A,Minv_array3(Q0));
    temp      = Mprod_array3(Mprod_array3(Minv_array3(Q0),Mpinv_array3(Astar)),b);
    
    % Contribution to ground reaction forces and moments and accelerations
    Contribution.muscleForce.(list{m+1}).F1R    = temp(1:3,:,:);
    Contribution.muscleForce.(list{m+1}).M1R    = Mprod_array3(Bstar1,temp(4:6,:,:));
    Contribution.muscleForce.(list{m+1}).d2Qdt2 = Mprod_array3(Mpinv_array3(Model.G),...
        Mprod_array3(Model.LG,...
        temp(7:30,:,:)));
end

% Compute contributions to the centre of pressure
% -------------------------------------------------------------------------
Y = [zeros(1,1,sw-1); ones(1,1,sw-1); zeros(1,1,sw-1)];
for i = 1:length(list)
    if i == 1
        F = Contribution.weight.F1R(:,:,2:sw);
        M = Contribution.weight.M1R(:,:,2:sw);
%         F(find(abs(F)<1)) = 0;
%         M(find(abs(M)<0.1)) = 0;
    else
        F = Contribution.muscleForce.(list{i}).F1R(:,:,2:sw);
        M = Contribution.muscleForce.(list{i}).M1R(:,:,2:sw);
%         F(find(abs(F)<1)) = 0;
%         M(find(abs(M)<0.1)) = 0;
    end
    temp = Mprod_array3(1/dot(F,F),cross(F,M)) - ...
           Mprod_array3(1/(dot(F,Y).*dot(F,F)), Mprod_array3(dot(F,M),cross(F,Y))) + ...
           Segment(1).Q(4:6,:,2:sw);
    lambda = Mprod_array3(-temp(2,:,:),1/F(2,:,:));
    CoP = Mprod_array3(1/dot(F,F),cross(F,M)) - ...
          Mprod_array3(1/(dot(F,Y).*dot(F,F)), Mprod_array3(dot(F,M),cross(F,Y))) + ...
          Segment(1).Q(4:6,:,2:sw) + Mprod_array3(lambda,F);
    CoPx(:,i) = permute(CoP(1,:,:),[3,2,1]);
    CoPz(:,i) = permute(CoP(3,:,:),[3,2,1]);
end

% CoP
% Proximal enpoint P (centre of pressure)
CoPx(:,12) = permute(Segment(1).Q(4,:,2:sw),[3,2,1]);
CoPz(:,12) = permute(Segment(1).Q(6,:,2:sw),[3,2,1]);
CoPx(isinf(CoPx)) = NaN;
CoPz(isinf(CoPz)) = NaN;

figure; hold on;
plot(CoPz(:,12),CoPx(:,12),'black')

%%
figure('units','normalized','position',[.15 .15 .3 .65]);
axis equal;
hold on;
% xlim([-0.13 0.31]);
% ylim([-0.03 0.39]);
plot(CoPz(:,1), CoPx(:,1),'Marker','.','Markersize',40,'Linestyle','none','Color',[0 113 188]/256);
plot(CoPz(:,i+1), CoPx(:,i+1),'Linestyle','-','Color','black','Linewidth',5);
% xlabel('Z_0 axis (m)');
ylabel('X_0 axis (m)');
set(gca,'FontSize',24,'FontName','Times New Roman');
set(findall(gcf,'type','text'),'FontSize',24,'FontName','Times New Roman');
cd('C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Contribution to CoP\Draft');
% saveas(gcf,'weight.pdf');

figure('units','normalized','position',[.15 .15 .3 .65]);
axis equal;
hold on;
% xlim([-0.13 0.31]);
% ylim([-0.03 0.39]);
plot(CoPz(:,2), CoPx(:,2),'Marker','.','Markersize',40,'Linestyle','none','Color',[216 82 24]/256);
plot(CoPz(:,3), CoPx(:,3),'Marker','.','Markersize',40,'Linestyle','none','Color',[236 176 31]/256);
plot(CoPz(:,i+1), CoPx(:,i+1),'Linestyle','-','Color','black','Linewidth',5);
% xlabel('Z_0 axis (m)');
% ylabel('X_0 axis (m)');
set(gca,'FontSize',24,'FontName','Times New Roman');
set(findall(gcf,'type','text'),'FontSize',24,'FontName','Times New Roman');
cd('C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Contribution to CoP\Draft');
% saveas(gcf,'hipAA.pdf');

figure('units','normalized','position',[.15 .15 .3 .65]);
axis equal;
hold on;
% xlim([-0.13 0.31]);
% ylim([-0.03 0.39]);
plot(CoPz(:,4), CoPx(:,4),'Marker','.','Markersize',40,'Linestyle','none','Color',[125 46 141]/256);
plot(CoPz(:,5), CoPx(:,5),'Marker','.','Markersize',40,'Linestyle','none','Color',[118 171 47]/256);
plot(CoPz(:,i+1), CoPx(:,i+1),'Linestyle','-','Color','black','Linewidth',5);
% xlabel('Z_0 axis (m)');
% ylabel('X_0 axis (m)');
set(gca,'FontSize',24,'FontName','Times New Roman');
set(findall(gcf,'type','text'),'FontSize',24,'FontName','Times New Roman');
cd('C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Contribution to CoP\Draft');
% saveas(gcf,'hipFE.pdf');

figure('units','normalized','position',[.15 .15 .3 .65]);
axis equal;
hold on;
% xlim([-0.13 0.31]);
% ylim([-0.03 0.39]);
plot(CoPz(:,6), CoPx(:,6),'Marker','.','Markersize',40,'Linestyle','none','Color',[24 56 230]/256);
plot(CoPz(:,7), CoPx(:,7),'Marker','.','Markersize',40,'Linestyle','none','Color',[100 176 240]/256);
plot(CoPz(:,i+1), CoPx(:,i+1),'Linestyle','-','Color','black','Linewidth',5);
xlabel('Z_0 axis (m)');
ylabel('X_0 axis (m)');
set(gca,'FontSize',24,'FontName','Times New Roman');
set(findall(gcf,'type','text'),'FontSize',24,'FontName','Times New Roman');
cd('C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Contribution to CoP\Draft');
% saveas(gcf,'kneeFE.pdf');

figure('units','normalized','position',[.15 .15 .3 .65]);
axis equal;
hold on;
% xlim([-0.13 0.31]);
% ylim([-0.03 0.39]);
plot(CoPz(:,8), CoPx(:,8),'Marker','.','Markersize',40,'Linestyle','none','Color',[125 46 141]/256);
plot(CoPz(:,9), CoPx(:,9),'Marker','.','Markersize',40,'Linestyle','none','Color',[118 171 47]/256);
plot(CoPz(:,i+1), CoPx(:,i+1),'Linestyle','-','Color','black','Linewidth',5);
xlabel('Z_0 axis (m)');
% ylabel('X_0 axis (m)');
set(gca,'FontSize',24,'FontName','Times New Roman');
set(findall(gcf,'type','text'),'FontSize',24,'FontName','Times New Roman');
cd('C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Contribution to CoP\Draft');
% saveas(gcf,'ankleFE.pdf');

figure('units','normalized','position',[.15 .15 .3 .65]);
axis equal;
hold on;
% xlim([-0.13 0.31]);
% ylim([-0.03 0.39]);
plot(CoPz(:,10), CoPx(:,10),'Marker','.','Markersize',40,'Linestyle','none','Color',[216 82 24]/256);
plot(CoPz(:,11), CoPx(:,11),'Marker','.','Markersize',40,'Linestyle','none','Color',[236 176 31]/256);
plot(CoPz(:,i+1), CoPx(:,i+1),'Linestyle','-','Color','black','Linewidth',5);
xlabel('Z_0 axis (m)');
% ylabel('X_0 axis (m)');
set(gca,'FontSize',24,'FontName','Times New Roman');
set(findall(gcf,'type','text'),'FontSize',24,'FontName','Times New Roman');
cd('C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Contribution to CoP\Draft');
% saveas(gcf,'ankleAA.pdf');

