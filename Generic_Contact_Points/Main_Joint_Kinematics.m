% MAIN PROGRAM
% Main_Joint_Kinematics.m
%__________________________________________________________________________
%
% PURPOSE
% Computation and plotting of 3D joint angles and displacements
%
% SYNOPSIS
% N/A (i.e., main program)
%
% DESCRIPTION
% Data loading, call of functions and plotting of joint coordinate system 
% angles and displacements
%
% REFERENCE
% R Dumas, T Robert, V Pomero, L Cheze. Joint and segment coordinate 
% systems revisited. Computer Methods in Biomechanics and Biomedical  
% Engineering. 2012;15(Suppl 1):183-5
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX)
% Mprod_array3.m
% Tinv_array3.m
% Q2Twu_array3.m
% Q2Tuv_array3.m
% Q2Tuw_array3.m
% R2mobileZXY_array3.m
% R2mobileZYX_array3
% Vnop_array3
% 
% MATLAB VERSION
% Matlab R2020a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas
% March 2010
%
% Modified by Raphaël Dumas
% October 2010
% Sequence ZYX for both ankle and wrist joints
% Figure captions
%
% Modified by Raphaël Dumas
% April 2012
% Version with 5 segments, i.e., including patella
%
% Modified by Raphael Dumas
% September 2012
% Normalisation of the 2d vector of the non-orthogonal projection on JCS axes
%
% Modified by Raphael Dumas
% March 2020
% Generic vs. Informed structures (Segment, Joint, Model)
% n,f,fc as Model.Informed fields
%__________________________________________________________________________


% Number of frames
n = Model.Informed.n;

% Interpolation parameters
k = (1:n)';
ko = (linspace(1,n,100))';

% Joint angles and displacements
for i = 2:5 % From i = 2 ankle to i = 5 hip
    
    % Transformation from the proximal segment axes
    % (with origin at endpoint D and with Z = w)
    % to the distal segment axes
    % (with origin at point P and with X = u)
    if i == 3 % Tibio-femoral
        Joint(3).Informed.T = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(5).Informed.Q)),... 
            Q2Tuv_array3(Segment(3).Informed.Q));
    elseif i == 5 % Hip
        Joint(5).Informed.T = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(6).Informed.Q)),...
            Q2Tuv_array3(Segment(5).Informed.Q));
        % Origin of proximal segment moved at the mean position of Pi
        % in proximal segment (rather than endpoint Di+1)
        Joint(5).Informed.T(1:3,4,:) = Joint(5).Informed.T(1:3,4,:) - ...
            repmat(mean(Joint(5).Informed.T(1:3,4,:),3),[1 1 n]);
    else % Ankle, patello-femoral
        Joint(i).Informed.T = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(i+1).Informed.Q)),...
            Q2Tuw_array3(Segment(i).Informed.Q));
    end
    
    % Joint coordinate system
    if i == 2 % ZYX sequence of mobile axis for ankle
        % Internal/extenal rotation on floating axis 
        % Euler angles
        Joint(i).Informed.Euler = R2mobileZYX_array3(Joint(i).Informed.T(1:3,1:3,:));
        % Joint displacement about the Euler angle axes
        Joint(i).Informed.dj = Vnop_array3(...
            Joint(i).Informed.T(1:3,4,:),... Di+1 to Pi in SCS of segment i+1
            repmat([0;0;1],[1 1 n]),... % Zi+1 in SCS of segment i+1
            Vnorm_array3(cross(repmat([0;0;1],[1 1 n]),Joint(i).Informed.T(1:3,1,:))),...
            Joint(i).Informed.T(1:3,1,:)); % Xi in SCS of segment i+1
    else % ZXY sequence of mobile axis
        % Euler angles
        Joint(i).Informed.Euler = R2mobileZXY_array3(Joint(i).Informed.T(1:3,1:3,:));
        % Joint displacement about the Euler angle axes
        Joint(i).Informed.dj = Vnop_array3(...
            Joint(i).Informed.T(1:3,4,:),... Di+1 to Pi in SCS of segment i+1
            repmat([0;0;1],[1 1 n]),... % Zi+1 in SCS of segment i+1
            Vnorm_array3(cross(Joint(i).Informed.T(1:3,2,:),repmat([0;0;1],[1 1 n]))),...
            Joint(i).Informed.T(1:3,2,:)); % Yi in SCS of segment i+1
    end

end

% 100% of gait cycle (or of propulsive cycle)
% Ankle joint angles and displacements
FE2 = interp1(k,permute(Joint(2).Informed.Euler(1,1,:),[3,2,1])*180/pi,ko,'spline');
IER2 = interp1(k,permute(Joint(2).Informed.Euler(1,2,:),[3,2,1])*180/pi,ko,'spline');
AA2 = interp1(k,permute(Joint(2).Informed.Euler(1,3,:),[3,2,1])*180/pi,ko,'spline');
LM2 = interp1(k,permute(Joint(2).Informed.dj(1,1,:),[3,2,1]),ko,'spline');
PD2 = interp1(k,permute(Joint(2).Informed.dj(2,1,:),[3,2,1]),ko,'spline');
AP2 = interp1(k,permute(Joint(2).Informed.dj(3,1,:),[3,2,1]),ko,'spline');
% Tibio-femoral joint angles and displacements
FE3 = interp1(k,permute(Joint(3).Informed.Euler(1,1,:),[3,2,1])*180/pi,ko,'spline');
AA3 = interp1(k,permute(Joint(3).Informed.Euler(1,2,:),[3,2,1])*180/pi,ko,'spline');
IER3 = interp1(k,permute(Joint(3).Informed.Euler(1,3,:),[3,2,1])*180/pi,ko,'spline');
LM3 = interp1(k,permute(Joint(3).Informed.dj(1,1,:),[3,2,1]),ko,'spline');
AP3 = interp1(k,permute(Joint(3).Informed.dj(2,1,:),[3,2,1]),ko,'spline');
PD3 = interp1(k,permute(Joint(3).Informed.dj(3,1,:),[3,2,1]),ko,'spline');
% Patello-femoral joint angles and displacements
FE4 = interp1(k,permute(Joint(4).Informed.Euler(1,1,:),[3,2,1])*180/pi,ko,'spline');
MLR4 = interp1(k,permute(Joint(4).Informed.Euler(1,2,:),[3,2,1])*180/pi,ko,'spline'); % Medial latetal rotation about X floating
MLT4 = interp1(k,permute(Joint(4).Informed.Euler(1,3,:),[3,2,1])*180/pi,ko,'spline'); % Medial lateral tilt about Y axis
LM4 = interp1(k,permute(Joint(4).Informed.dj(1,1,:),[3,2,1]),ko,'spline'); % Lateral medial shift
AP4 = interp1(k,permute(Joint(4).Informed.dj(2,1,:),[3,2,1]),ko,'spline');
PD4 = interp1(k,permute(Joint(4).Informed.dj(3,1,:),[3,2,1]),ko,'spline');
% Hip joint angles and displacements
FE5 = interp1(k,permute(Joint(5).Informed.Euler(1,1,:),[3,2,1])*180/pi,ko,'spline');
AA5 = interp1(k,permute(Joint(5).Informed.Euler(1,2,:),[3,2,1])*180/pi,ko,'spline');
IER5 = interp1(k,permute(Joint(5).Informed.Euler(1,3,:),[3,2,1])*180/pi,ko,'spline');
LM5 = interp1(k,permute(Joint(5).Informed.dj(1,1,:),[3,2,1]),ko,'spline');
AP5 = interp1(k,permute(Joint(5).Informed.dj(2,1,:),[3,2,1]),ko,'spline');
PD5 = interp1(k,permute(Joint(5).Informed.dj(3,1,:),[3,2,1]),ko,'spline');


% Figure for ankle
figure;
hold on;
% Flexion Extension
subplot(2,3,1);
hold on;
plot(FE2,'LineWidth',1,'Color','k');
title ('Right Ankle Flexion (+) / Extension (-)');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Adduction Abduction
subplot(2,3,2);
hold on;
plot(AA2,'LineWidth',1,'Color','k');
title ('Right Ankle Adduction (+) / Abduction (-)');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Internal External Rotation
subplot(2,3,3);
hold on;
plot(IER2,'LineWidth',1,'Color','k');
title ('Right Ankle Internal (+) / External (-) Rotation');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Lateral Medial
subplot(2,3,4);
hold on;
plot(LM2*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Ankle Lateral (+) /  Medial (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');
% Anterior Posterior
subplot(2,3,5);
hold on;
plot(AP2*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Ankle Anterior (+) / Posterior (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');
% Proximal Distal
subplot(2,3,6);
hold on;
plot(PD2*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Ankle Proximal (+) / Distal (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');


% Figures for tibio-femoral
figure;
hold on;
% Extension Flexion
subplot(2,3,1);
hold on;
plot(FE3,'LineWidth',1,'Color','k');
title ('Right Tibio-Femoral Extension (+) / Flexion (-)');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Adduction Abduction
subplot(2,3,2);
hold on;
plot(AA3,'LineWidth',1,'Color','k');
title ('Right Tibio-Femoral Adduction (+) / Abduction (-)');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Internal External Rotation
subplot(2,3,3);
hold on;
plot(IER3,'LineWidth',1,'Color','k');
title ('Right Tibio-Femoral Internal (+) / External (-) Rotation');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Lateral Medial
subplot(2,3,4);
hold on;
plot(LM3*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Tibio-Femoral Lateral (+) /  Medial (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');
% Anterior Posterior
subplot(2,3,5);
hold on;
plot(AP3*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Tibio-Femoral Anterior (+) / Posterior (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');
% Proximal Distal
subplot(2,3,6);
hold on;
plot(PD3*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Tibio-Femoral Proximal (+) / Distal (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');


% Figures for patello-femoral
figure;
hold on;
% Extension Flexion
subplot(2,3,1);
hold on;
plot(FE4,'LineWidth',1,'Color','k');
title ('Right Patello-Femoral Extension (+) / Flexion (-)');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Medial latetal rotation
subplot(2,3,2);
hold on;
plot(MLR4,'LineWidth',1,'Color','k');
title ('Right Patello-Femoral Medial (+) / Latetal (-) Rotation');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Medial latetal tilt
subplot(2,3,3);
hold on;
plot(MLT4,'LineWidth',1,'Color','k');
title ('Right Patello-Femoral Medial (+) / Latetal (-) Tilt');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Lateral Medial
subplot(2,3,4);
hold on;
plot(LM4*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Patello-Femoral Lateral (+) /  Medial (-) Shift');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');
% Anterior Posterior
subplot(2,3,5);
hold on;
plot(AP4*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Patello-Femoral Anterior (+) / Posterior (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');
% Proximal Distal
subplot(2,3,6);
hold on;
plot(PD4*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Patello-Femoral Proximal (+) / Distal (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');

% Figure for hip
figure;
hold on;
% Flexion Extension 
subplot(2,3,1);
hold on;
plot(FE5,'LineWidth',1,'Color','k');
title ('Right Hip Flexion (+) / Extension (-)');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Adduction Abduction
subplot(2,3,2);
hold on;
plot(AA5,'LineWidth',1,'Color','k');
title ('Right Hip Adduction (+) / Abduction (-)');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Internal External Rotation
subplot(2,3,3);
hold on;
plot(IER5,'LineWidth',1,'Color','k');
title ('Right Hip Internal (+) / External (-) Rotation');
xlabel('% of Gait Cycle');
ylabel('Angle (in degree)');
% Lateral Medial
subplot(2,3,4);
hold on;
plot(LM5*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Hip Lateral (+) /  Medial (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');
% Anterior Posterior
subplot(2,3,5);
hold on;
plot(AP5*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Hip Anterior (+) / Posterior (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');
% Proximal Distal
subplot(2,3,6);
hold on;
plot(PD5*1000,'LineWidth',1,'Color','k'); % mm
title ('Right Hip Proximal (+) / Distal (-)');
xlabel('% of Gait Cycle');
ylabel('Displacement (in mm)');
