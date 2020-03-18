% MAIN PROGRAM
% Main_Muscle_Lines_Visualisation_PS.m
%__________________________________________________________________________
%
% PURPOSE
% Visualisation of muscle lines of action
%
% SYNOPSIS
% N/A (i.e., main program)
%
% DESCRIPTION
% Plotting of segment and muscle lines of actions
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
%  Mprod_array3.m
%
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphael Dumas
% April 2018
%
% Modified by Raphael Dumas
% June 2019
% Visualisation of prosthetic component
% Plot of contact points
% Renamed _GC5
%
% Modified by Raphael Dumas
% January 2020
% Renamed _PS
%__________________________________________________________________________

% Number of frames
n = size(Segment(2).Q,3);

% Frames of interest
k = round(55*length(Segment(2).Q)/100);

figure
hold on
axis equal

% -------------------------------------------------------------------------
% Stl files
% CAD from PS
% -------------------------------------------------------------------------

% Virtual palpation in CT-scan/CAD model (tibia)
Tt_implant2bone = [ 0.9999    0.0112   -0.0093   -0.0027; ...
   -0.0111    0.9999    0.0105   -0.0043; ...
    0.0094   -0.0104    0.9999   -0.0027; ...
         0         0         0    1.0000];
     
% From CT-scan/CAD femur origin to CT-scan/CAD tibia origin
Tfimplant2timplant = [eye(3),-[0; 43.06; 0]/1000; [0,0,0,1]];

% Virtual palpation in CT-scan/CAD model (femur)
Tf_implant2bone = [0.9999    0.0048   -0.0096   -0.0027; ...
   -0.0049    0.9999   -0.0135   -0.0043; ...
    0.0095    0.0135    0.9999   -0.0027; ...
         0         0         0    1.0000]; 

for i = 1:3
    switch i
        case 1
            [Triangles,Points] = stlread('Tibial_Insert.stl');
            C = 'b';
            Points(:,3) = - Points(:,3); % Left knee
            Segment(3).Vertices = (inv(Tt_implant2bone)*Tfimplant2timplant*...
                [Points'/1000;ones(1,length(Points))])';
            Segment(3).Vertices(:,4) = [];
            Segment(3).Faces = Triangles;        
        case 2
            [Triangles,Points] = stlread('Tibia.stl');
            C = 'k';
            Points(:,3) = - Points(:,3); % Left knee
            Segment(3).Vertices = (inv(Tt_implant2bone)*Tfimplant2timplant*...
                [Points'/1000;ones(1,length(Points))])';
            Segment(3).Vertices(:,4) = [];
            Segment(3).Faces = Triangles;        
        case 3
            [Triangles,Points] = stlread('Fibula.stl');
            C = 'k';
            Points(:,3) = - Points(:,3); % Left knee
            Segment(3).Vertices = (inv(Tt_implant2bone)*Tfimplant2timplant*...
                [Points'/1000;ones(1,length(Points))])';
            Segment(3).Vertices(:,4) = [];
            Segment(3).Faces = Triangles;        
    end
    Segment(3).T = Q2Tuv_array3(Segment(3).Q);
    Vertices_k = (Segment(3).T(:,:,k)*[Segment(3).Vertices, ...
        ones(length(Segment(3).Vertices),1)]')';
    patch('Faces',Segment(3).Faces,'Vertices',Vertices_k(:,1:3),...
        'FaceColor',C,'FaceAlpha', 0.15,'EdgeColor','none');
end

for i = 1:2
    switch i
        case 1
            [Triangles,Points] = stlread('Femoral_Component.stl');
            C = 'b';
            Points(:,3) = - Points(:,3); % Left knee
            Segment(5).Vertices = [0, - Segment(5).L, 0, 0] + ... Origin was at rD5
                (inv(Tf_implant2bone)*[Points'/1000;ones(1,length(Points))])';
            Segment(5).Vertices(:,4) = [];
            Segment(5).Faces = Triangles;
        case 2
            [Triangles,Points] = stlread('Femur.stl');
            C = 'k';
            Points(:,3) = - Points(:,3); % Left knee
            Segment(5).Vertices = [0, - Segment(5).L, 0, 0] + ... Origin was at rD5
                (inv(Tf_implant2bone)*[Points'/1000;ones(1,length(Points))])';
            Segment(5).Vertices(:,4) = [];
            Segment(5).Faces = Triangles;
    end
    Segment(5).T = Q2Tuv_array3(Segment(5).Q);
    Vertices_k = (Segment(5).T(:,:,k)*[Segment(5).Vertices, ...
        ones(length(Segment(5).Vertices),1)]')';
    patch('Faces',Segment(5).Faces,'Vertices',Vertices_k(:,1:3),...
        'FaceColor',C,'FaceAlpha', 0.15,'EdgeColor','none');
end


% -------------------------------------------------------------------------
% Interpolation matrices
% Insertion and origin of muscle lines of action
% -------------------------------------------------------------------------

for i = 2:6
    switch i
        case 2 % Virtual markers (3 first are for ankle joint)
            nj = 4:15;
        case 3 % Virtual markers (9 first are for ankle (), tibio-femoral and patello-femoral joints)
            nj = 9:28; % V93 is patettar tendon insertion
        case 4
            nj = 2:5; % Virtual markers (first is for patello-femoral joint)
        case 5
            nj = 6:43; % Virtual markers (5 first are for tibio-femoral joint)
        case 6
            nj = 2:28; % Virtual markers (first is for hip joint)
    end
    for j = nj
        % Interpolation matrix
        eval(['NV',num2str(j),num2str(i),...
            '= [Segment(',num2str(i),').nV(1,',num2str(j),')*eye(3),',...
            '(1 + Segment(',num2str(i),').nV(2,',num2str(j),'))*eye(3),',...
            '- Segment(',num2str(i),').nV(2,',num2str(j),')*eye(3),',...
            'Segment(',num2str(i),').nV(3,',num2str(j),')*eye(3)];']);
        % Point of insertion and origin lines of muscle action
        eval(['Segment(',num2str(i),').rV',num2str(j),num2str(i),...
            '(1:3,1,:) = Mprod_array3(repmat(NV',num2str(j),num2str(i),...
            ',[1,1,n]),Segment(',num2str(i),').Q);']);
    end
end

% -------------------------------------------------------------------------
% Segment and lines of action
% -------------------------------------------------------------------------

% PLOT SEGMENT
% Segment 2 'Foot'
plot3([Segment(2).Q(4,1,k),Segment(2).Q(7,1,k)],[Segment(2).Q(5,1,k),Segment(2).Q(8,1,k)],[Segment(2).Q(6,1,k),Segment(2).Q(9,1,k)],'k','linewidth',3);
% Segment 3 'Shank'
plot3([Segment(3).Q(4,1,k),Segment(3).Q(7,1,k)],[Segment(3).Q(5,1,k),Segment(3).Q(8,1,k)],[Segment(3).Q(6,1,k),Segment(3).Q(9,1,k)],'k','linewidth',3);
% Segment 4 'Patella'
plot3([Segment(4).Q(4,1,k),Segment(4).Q(7,1,k)],[Segment(4).Q(5,1,k),Segment(4).Q(8,1,k)],[Segment(4).Q(6,1,k),Segment(4).Q(9,1,k)],'k','linewidth',3);
% Segment 5 'Femur'
plot3([Segment(5).Q(4,1,k),Segment(5).Q(7,1,k)],[Segment(5).Q(5,1,k),Segment(5).Q(8,1,k)],[Segment(5).Q(6,1,k),Segment(5).Q(9,1,k)],'k','linewidth',3);
% Segment 6 'Pelvis'
plot3([Segment(6).Q(4,1,k),Segment(6).Q(7,1,k)],[Segment(6).Q(5,1,k),Segment(6).Q(8,1,k)],[Segment(6).Q(6,1,k),Segment(6).Q(9,1,k)],'k','linewidth',3);

% PLOT MUSCLE
% Gluteus maximus I
plot3([Segment(6).rV26(1,1,k),Segment(5).rV75(1,1,k)],[Segment(6).rV26(2,1,k),Segment(5).rV75(2,1,k)],[Segment(6).rV26(3,1,k),Segment(5).rV75(3,1,k)],'r','linewidth',2);
% Gluteus maximus II
plot3([Segment(6).rV36(1,1,k),Segment(5).rV85(1,1,k)],[Segment(6).rV36(2,1,k),Segment(5).rV85(2,1,k)],[Segment(6).rV36(3,1,k),Segment(5).rV85(3,1,k)],'r','linewidth',2);
% Gluteus maximus III
plot3([Segment(6).rV46(1,1,k),Segment(5).rV95(1,1,k)],[Segment(6).rV46(2,1,k),Segment(5).rV95(2,1,k)],[Segment(6).rV46(3,1,k),Segment(5).rV95(3,1,k)],'r','linewidth',2);
% Gluteus medius I
plot3([Segment(6).rV56(1,1,k),Segment(5).rV105(1,1,k)],[Segment(6).rV56(2,1,k),Segment(5).rV105(2,1,k)],[Segment(6).rV56(3,1,k),Segment(5).rV105(3,1,k)],'r','linewidth',2);
% Gluteus medius II
plot3([Segment(6).rV66(1,1,k),Segment(5).rV115(1,1,k)],[Segment(6).rV66(2,1,k),Segment(5).rV115(2,1,k)],[Segment(6).rV66(3,1,k),Segment(5).rV115(3,1,k)],'r','linewidth',2);
% Gluteus medius III
plot3([Segment(6).rV76(1,1,k),Segment(5).rV125(1,1,k)],[Segment(6).rV76(2,1,k),Segment(5).rV125(2,1,k)],[Segment(6).rV76(3,1,k),Segment(5).rV125(3,1,k)],'r','linewidth',2);
% Gluteus minimus I
plot3([Segment(6).rV86(1,1,k),Segment(5).rV135(1,1,k)],[Segment(6).rV86(2,1,k),Segment(5).rV135(2,1,k)],[Segment(6).rV86(3,1,k),Segment(5).rV135(3,1,k)],'r','linewidth',2);
% Gluteus minimus II
plot3([Segment(6).rV96(1,1,k),Segment(5).rV145(1,1,k)],[Segment(6).rV96(2,1,k),Segment(5).rV145(2,1,k)],[Segment(6).rV96(3,1,k),Segment(5).rV145(3,1,k)],'r','linewidth',2);
% Gluteus minimus III
plot3([Segment(6).rV106(1,1,k),Segment(5).rV155(1,1,k)],[Segment(6).rV106(2,1,k),Segment(5).rV155(2,1,k)],[Segment(6).rV106(3,1,k),Segment(5).rV155(3,1,k)],'r','linewidth',2);
% Adductor longus
plot3([Segment(6).rV116(1,1,k),Segment(5).rV165(1,1,k)],[Segment(6).rV116(2,1,k),Segment(5).rV165(2,1,k)],[Segment(6).rV116(3,1,k),Segment(5).rV165(3,1,k)],'r','linewidth',2);
% Adductor brevis
plot3([Segment(6).rV126(1,1,k),Segment(5).rV175(1,1,k)],[Segment(6).rV126(2,1,k),Segment(5).rV175(2,1,k)],[Segment(6).rV126(3,1,k),Segment(5).rV175(3,1,k)],'r','linewidth',2);
% Adductor magnus I
plot3([Segment(6).rV136(1,1,k),Segment(5).rV185(1,1,k)],[Segment(6).rV136(2,1,k),Segment(5).rV185(2,1,k)],[Segment(6).rV136(3,1,k),Segment(5).rV185(3,1,k)],'r','linewidth',2);
% Adductor magnus II
plot3([Segment(6).rV146(1,1,k),Segment(5).rV195(1,1,k)],[Segment(6).rV146(2,1,k),Segment(5).rV195(2,1,k)],[Segment(6).rV146(3,1,k),Segment(5).rV195(3,1,k)],'r','linewidth',2);
% Adductor magnus III
plot3([Segment(6).rV156(1,1,k),Segment(5).rV205(1,1,k)],[Segment(6).rV156(2,1,k),Segment(5).rV205(2,1,k)],[Segment(6).rV156(3,1,k),Segment(5).rV205(3,1,k)],'r','linewidth',2);
% Pectineus
plot3([Segment(6).rV166(1,1,k),Segment(5).rV215(1,1,k)],[Segment(6).rV166(2,1,k),Segment(5).rV215(2,1,k)],[Segment(6).rV166(3,1,k),Segment(5).rV215(3,1,k)],'r','linewidth',2);
% Illiacus
plot3([Segment(6).rV176(1,1,k),Segment(5).rV225(1,1,k)],[Segment(6).rV176(2,1,k),Segment(5).rV225(2,1,k)],[Segment(6).rV176(3,1,k),Segment(5).rV225(3,1,k)],'r','linewidth',2);
% Psoas
plot3([Segment(6).rV186(1,1,k),Segment(5).rV235(1,1,k)],[Segment(6).rV186(2,1,k),Segment(5).rV235(2,1,k)],[Segment(6).rV186(3,1,k),Segment(5).rV235(3,1,k)],'r','linewidth',2);
% Quadratus femoris
plot3([Segment(6).rV196(1,1,k),Segment(5).rV245(1,1,k)],[Segment(6).rV196(2,1,k),Segment(5).rV245(2,1,k)],[Segment(6).rV196(3,1,k),Segment(5).rV245(3,1,k)],'r','linewidth',2);
% Gemellus
plot3([Segment(6).rV206(1,1,k),Segment(5).rV255(1,1,k)],[Segment(6).rV206(2,1,k),Segment(5).rV255(2,1,k)],[Segment(6).rV206(3,1,k),Segment(5).rV255(3,1,k)],'r','linewidth',2);
% Piriformis
plot3([Segment(6).rV216(1,1,k),Segment(5).rV265(1,1,k)],[Segment(6).rV216(2,1,k),Segment(5).rV265(2,1,k)],[Segment(6).rV216(3,1,k),Segment(5).rV265(3,1,k)],'r','linewidth',2);
% Tensor fasciae latae
plot3([Segment(6).rV226(1,1,k),Segment(5).rV275(1,1,k),Segment(5).rV285(1,1,k),Segment(3).rV103(1,1,k)],...
    [Segment(6).rV226(2,1,k),Segment(5).rV275(2,1,k),Segment(5).rV285(2,1,k),Segment(3).rV103(2,1,k)],...
    [Segment(6).rV226(3,1,k),Segment(5).rV275(3,1,k),Segment(5).rV285(3,1,k),Segment(3).rV103(3,1,k)],'r','linewidth',2);
% Gracilis
plot3([Segment(6).rV236(1,1,k),Segment(3).rV113(1,1,k)],[Segment(6).rV236(2,1,k),Segment(3).rV113(2,1,k)],[Segment(6).rV236(3,1,k),Segment(3).rV113(3,1,k)],'r','linewidth',2);
% Sartorius
plot3([Segment(6).rV276(1,1,k),Segment(5).rV295(1,1,k),Segment(3).rV153(1,1,k)],...
    [Segment(6).rV276(2,1,k),Segment(5).rV295(2,1,k),Segment(3).rV153(2,1,k)],...
    [Segment(6).rV276(3,1,k),Segment(5).rV295(3,1,k),Segment(3).rV153(3,1,k)],'r','linewidth',2);
% Semimenbranosus
plot3([Segment(6).rV246(1,1,k),Segment(3).rV123(1,1,k)],[Segment(6).rV246(2,1,k),Segment(3).rV123(2,1,k)],[Segment(6).rV246(3,1,k),Segment(3).rV123(3,1,k)],'r','linewidth',2);
% Semitendinosus
plot3([Segment(6).rV256(1,1,k),Segment(3).rV133(1,1,k)],[Segment(6).rV256(2,1,k),Segment(3).rV133(2,1,k)],[Segment(6).rV256(3,1,k),Segment(3).rV133(3,1,k)],'r','linewidth',2);
% Biceps femoris long head
plot3([Segment(6).rV266(1,1,k),Segment(3).rV143(1,1,k)],[Segment(6).rV266(2,1,k),Segment(3).rV143(2,1,k)],[Segment(6).rV266(3,1,k),Segment(3).rV143(3,1,k)],'r','linewidth',2);
% Biceps femoris short head
plot3([Segment(5).rV305(1,1,k),Segment(3).rV163(1,1,k)],[Segment(5).rV305(2,1,k),Segment(3).rV163(2,1,k)],[Segment(5).rV305(3,1,k),Segment(3).rV163(3,1,k)],'r','linewidth',2);
% Rectus femoris (in case of tibio-femoral angle > - 83.65 degrees)
plot3([Segment(6).rV286(1,1,k),Segment(4).rV24(1,1,k)],[Segment(6).rV286(2,1,k),Segment(4).rV24(2,1,k)],[Segment(6).rV286(3,1,k),Segment(4).rV24(3,1,k)],'r','linewidth',2);
% Vastus medialis (in case of tibio-femoral angle > -69.32 degrees)
plot3([Segment(5).rV325(1,1,k),Segment(4).rV34(1,1,k)],[Segment(5).rV325(2,1,k),Segment(4).rV34(2,1,k)],[Segment(5).rV325(3,1,k),Segment(4).rV34(3,1,k)],'r','linewidth',2);
% Vastus intermedialis (in case of tibio-femoral angle > - -81.36 degrees)
plot3([Segment(5).rV355(1,1,k),Segment(4).rV44(1,1,k)],[Segment(5).rV355(2,1,k),Segment(4).rV44(2,1,k)],[Segment(5).rV355(3,1,k),Segment(4).rV44(3,1,k)],'r','linewidth',2);
% Vastus lateralis (in case of tibio-femoral angle > -69.32 degrees)
plot3([Segment(5).rV375(1,1,k),Segment(4).rV54(1,1,k)],[Segment(5).rV375(2,1,k),Segment(4).rV54(2,1,k)],[Segment(5).rV375(3,1,k),Segment(4).rV54(3,1,k)],'r','linewidth',2);
% Gastrocnemius medialis (in case of tibio-femoral angle > -44.12 degrees)
plot3([Segment(5).rV405(1,1,k),Segment(3).rV173(1,1,k),Segment(2).rV42(1,1,k)],...
    [Segment(5).rV405(2,1,k),Segment(3).rV173(2,1,k),Segment(2).rV42(2,1,k)],...
    [Segment(5).rV405(3,1,k),Segment(3).rV173(3,1,k),Segment(2).rV42(3,1,k)],'r','linewidth',2);
% Gastrocnemius lateralis (in case of tibio-femoral angle > -44.12 degrees)
plot3([Segment(5).rV425(1,1,k),Segment(3).rV183(1,1,k),Segment(2).rV52(1,1,k)],...
    [Segment(5).rV425(2,1,k),Segment(3).rV183(2,1,k),Segment(2).rV52(2,1,k)],...
    [Segment(5).rV425(3,1,k),Segment(3).rV183(3,1,k),Segment(2).rV52(3,1,k)],'r','linewidth',2);
% Soleus
plot3([Segment(3).rV193(1,1,k),Segment(2).rV62(1,1,k)],[Segment(3).rV193(2,1,k),Segment(2).rV62(2,1,k)],[Segment(3).rV193(3,1,k),Segment(2).rV62(3,1,k)],'r','linewidth',2);
% Tibialis posterior
plot3([Segment(3).rV203(1,1,k),Segment(2).rV72(1,1,k)],[Segment(3).rV203(2,1,k),Segment(2).rV72(2,1,k)],[Segment(3).rV203(3,1,k),Segment(2).rV72(3,1,k)],'r','linewidth',2);
% Tibialis anterior
plot3([Segment(3).rV213(1,1,k),Segment(2).rV82(1,1,k)],[Segment(3).rV213(2,1,k),Segment(2).rV82(2,1,k)],[Segment(3).rV213(3,1,k),Segment(2).rV82(3,1,k)],'r','linewidth',2);
% Peroneus brevis
plot3([Segment(3).rV223(1,1,k),Segment(2).rV92(1,1,k)],[Segment(3).rV223(2,1,k),Segment(2).rV92(2,1,k)],[Segment(3).rV223(3,1,k),Segment(2).rV92(3,1,k)],'r','linewidth',2);
% Peroneus longus
plot3([Segment(3).rV233(1,1,k),Segment(2).rV102(1,1,k)],[Segment(3).rV233(2,1,k),Segment(2).rV102(2,1,k)],[Segment(3).rV233(3,1,k),Segment(2).rV102(3,1,k)],'r','linewidth',2);
% Peroneus tertius
plot3([Segment(3).rV243(1,1,k),Segment(2).rV112(1,1,k)],[Segment(3).rV243(2,1,k),Segment(2).rV112(2,1,k)],[Segment(3).rV243(3,1,k),Segment(2).rV112(3,1,k)],'r','linewidth',2);
% Extensor digitorum longus
plot3([Segment(3).rV253(1,1,k),Segment(2).rV122(1,1,k)],[Segment(3).rV253(2,1,k),Segment(2).rV122(2,1,k)],[Segment(3).rV253(3,1,k),Segment(2).rV122(3,1,k)],'r','linewidth',2);
% Extensor hallucis longus
plot3([Segment(3).rV263(1,1,k),Segment(2).rV132(1,1,k)],[Segment(3).rV263(2,1,k),Segment(2).rV132(2,1,k)],[Segment(3).rV263(3,1,k),Segment(2).rV132(3,1,k)],'r','linewidth',2);
% Flexor digitorum longus
plot3([Segment(3).rV273(1,1,k),Segment(2).rV142(1,1,k)],[Segment(3).rV273(2,1,k),Segment(2).rV142(2,1,k)],[Segment(3).rV273(3,1,k),Segment(2).rV142(3,1,k)],'r','linewidth',2);
% Flexor hallucis longus
plot3([Segment(3).rV283(1,1,k),Segment(2).rV152(1,1,k)],[Segment(3).rV283(2,1,k),Segment(2).rV152(2,1,k)],[Segment(3).rV283(3,1,k),Segment(2).rV152(3,1,k)],'r','linewidth',2);

legend({'Tibial Insert', 'Tibia', 'Fibula', ...
    'Femoral Component', 'Femur', ...
    'Foot','Shank','Patella','Thigh','Pelvis', ...
    'Gluteus maximus I', 'Gluteus maximus II', 'Gluteus maximus III', ...
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

% -------------------------------------------------------------------------
% More plots
% -------------------------------------------------------------------------

% PLOT JOINT
% % Contact points
plot3(Segment(3).rV43(1,1,k),Segment(3).rV43(2,1,k),Segment(3).rV43(3,1,k),'ob');
plot3(Segment(3).rV53(1,1,k),Segment(3).rV53(2,1,k),Segment(3).rV53(3,1,k),'ob');
% Patellofemoral tendon
plot3([Segment(3).rV93(1,1,k),Segment(4).Q(7,1,k)],[Segment(3).rV93(2,1,k),Segment(4).Q(8,1,k)],[Segment(3).rV93(3,1,k),Segment(4).Q(9,1,k)],'b','linewidth',2);
% Patellofemoral axis point
plot3(Segment(5).rV65(1,1,k),Segment(5).rV65(2,1,k),Segment(5).rV65(3,1,k),'ob');

% Markers
for i = [2,3,5,6]
    for j = 1:size(Segment(i).rM,2)
        plot3(permute(Segment(i).rM(1,j,k),[3,2,1]),...
            permute(Segment(i).rM(2,j,k),[3,2,1]),...
            permute(Segment(i).rM(3,j,k),[3,2,1]),'ok');
    end
end

