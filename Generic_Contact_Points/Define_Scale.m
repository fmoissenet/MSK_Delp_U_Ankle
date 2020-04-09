%  FUNCTION
% Define_Scale.m
%__________________________________________________________________________
%
% PURPOSE
% Define scales to be applied to generic model data
%
% SYNOPSIS
% [Segment,Model] = Define_Scale(Segment,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model 
%
% OUTPUT
% Segment (cf. data structure in user guide)
% Model
%
% DESCRIPTION
% Define scales for all segments and for the whole model based on segment
% lengths with kept generic proportions 
%
% REFERENCE
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% 
% MATLAB VERSION
% Matlab R2020a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphael Dumas
% March 2020
%__________________________________________________________________________

function [Segment,Model] = Define_Scale(Segment,Model)

% Number of frames
n = Model.Informed.n;


%% -------------------------------------------------------------------------
% Mean segment geometry (by default)
% -------------------------------------------------------------------------

for i = [2,3,5,6] % From i = 2 (Foot) to i = 6 (Pelvis)
    
    % Segment length
    Segment(i).Informed.L = mean(sqrt(sum((Segment(i).Informed.Q(4:6,1,:) - ...
            Segment(i).Informed.Q(7:9,1,:)).^2)),3);
    
    % Alpha angle between (rP - rD) and w
    Segment(i).Informed.alpha = mean(acosd(dot(Segment(i).Informed.Q(4:6,1,:) - ...
        Segment(i).Informed.Q(7:9,1,:), Segment(i).Informed.Q(10:12,1,:))./...
        sqrt(sum((Segment(i).Informed.Q(4:6,1,:) - ...
        Segment(i).Informed.Q(7:9,1,:)).^2))),3);
    
    % Beta angle between u and w
    Segment(i).Informed.beta = mean(acosd(dot(Segment(i).Informed.Q(10:12,1,:), ...
        Segment(i).Informed.Q(1:3,1,:))),3);
    
    % Gamma angle between u and (rP - rD)
    Segment(i).Informed.gamma = mean(acosd(dot(Segment(i).Informed.Q(1:3,1,:), ...
        Segment(i).Informed.Q(4:6,1,:) - Segment(i).Informed.Q(7:9,1,:))./...
        sqrt(sum((Segment(i).Informed.Q(4:6,1,:) - ...
        Segment(i).Informed.Q(7:9,1,:)).^2))),3);
    
    % Transformation from ICS to SCS
    Segment(i).Informed.T = Q2Tuv_array3(Segment(i).Informed.Q);
    
end

%% -------------------------------------------------------------------------
% Built averaged and aligned informed model
% -------------------------------------------------------------------------

% Pelvis
Informed_L6 = Segment(6).Informed.L;
Informed_gamma6 = Segment(6).Informed.gamma;
Informed_rD6 = [0;0;0];
Informed_rP6 = Informed_rD6 + ...
    Informed_L6.*cosd(Informed_gamma6)*[1;0;0] + ... % About X
    Informed_L6.*sind(Informed_gamma6)*[0;1;0]; % About Y
Informed_nV16 = mean(Vnop_array3(...
    Segment(5).Informed.Q(4:6,1,:) - Segment(6).Informed.Q(4:6,1,:),...
    Segment(6).Informed.Q(1:3,1,:),...
    Segment(6).Informed.Q(4:6,1,:) - Segment(6).Informed.Q(7:9,1,:),...
    Segment(6).Informed.Q(10:12,1,:)),3);

% Thigh
Informed_L5 = Segment(5).Informed.L;
Informed_rP5 = Informed_rD6 + ...
    Informed_nV16(3,1)*[0;0;1]; % About Z
Informed_rD5 = Informed_rP5 - ...
    Informed_L5*[0;1;0]; % About Y

% Shank
Informed_L3 = Segment(3).Informed.L;
Informed_rD3 = Informed_rD5 - ... % rP3 = rD5
    Informed_L3*[0;1;0]; % About Y

% Foot
Informed_L2 = Segment(2).Informed.L;
Informed_gamma2 = Segment(2).Informed.gamma;
Informed_rD2 = Informed_rD3 - ... % rP2 = rD3
    Informed_L2.*cosd(Informed_gamma2)*[1;0;0] - ... % About X
    Informed_L2.*sind(Informed_gamma2)*[0;1;0]; % About Y


%% -------------------------------------------------------------------------
% Built aligned generic model
% -------------------------------------------------------------------------

% Pelvis
Generic_L6 = Segment(6).Generic.L;
Generic_gamma6 = Segment(6).Generic.gamma;
Generic_rD6 = [0;0;0];
Generic_rP6 = Informed_rD6 + ...
    Generic_L6.*cosd(Generic_gamma6)*[1;0;0] + ... % About X
    Generic_L6.*sind(Generic_gamma6)*[0;1;0]; % About Y
Generic_rV16 = Segment(6).Generic.rVs(:,1);

% Thigh
Generic_L5 = Segment(5).Generic.L;
Generic_rP5 = Generic_rD6 + ...
    Generic_rV16(3,1)*[0;0;1]; % About Z
Generic_rD5 = Generic_rP5 - ...
    Generic_L5*[0;1;0]; % About Y

% Shank
Generic_L3 = Segment(3).Generic.L;
Generic_rD3 = Generic_rD5 - ... % rP3 = rD5
    Generic_L3*[0;1;0]; % About Y

% Foot
Generic_L2 = Segment(2).Generic.L;
Generic_gamma2 = Segment(2).Generic.gamma;
Generic_rD2 = Generic_rD3 - ... % rP2 = rD3
    Generic_L2.*cosd(Generic_gamma2)*[1;0;0] - ... % About X
    Generic_L2.*sind(Generic_gamma2)*[0;1;0]; % About Y


%% -------------------------------------------------------------------------
% Best scale, anterior-posterior, and superior-inferior translations 
% -------------------------------------------------------------------------

% Generic control points
Xi = [Generic_rD6,Generic_rP6,...
    Generic_rP5,Generic_rD5,Generic_rD3,Generic_rD2, ...
    Generic_rD5,Generic_rD3]'; % Double these endpoints (most trustable)

% Informed control points
Xj = [Informed_rD6,Informed_rP6,...
    Informed_rP5,Informed_rD5,Informed_rD3,Informed_rD2, ...
    Informed_rD5,Informed_rD3]'; % Double these endpoints (most trustable)

% Minimum of sum of squared errors
C = @(x)SSE(x,Xi,Xj); % Anonymous function (to pass extra parameters)
x0 = [1;0;0];
xmin = [0.5,-0.005,-0.005]; % half scale, 5 mm
xmax = [2,0.005,0.005]; % double scale, 5 mm
x = fmincon(C,x0,[],[],[],[],xmin,xmax) % Display
A(2:4,1:3) = eye(3)*x(1,1); % Scale
A(1,1) = x(2,1); % Anterior-posterior translation
A(1,2) = x(3,1); % Superior-inferior translation
% Final control points minus transformed initial control points
Xjs =  Xj - [ones(size(Xi,1),1),Xi]*A; % Apply scale
% Sum of squared errors (residual)
diag(Xjs'*Xjs) % Display SSE about each axis

% Uniform scale
Model.Informed.Scale = x(1,1);

% Scaled generic aligned model
Scaled_rP6 = ([1;Generic_rP6]'*A)';
Scaled_rD6 = ([1;Generic_rD6]'*A)';
Scaled_rP5 = ([1;Generic_rP5]'*A)';
Scaled_rD5 = ([1;Generic_rD5]'*A)';
Scaled_rD3 = ([1;Generic_rD3]'*A)';
Scaled_rD2 = ([1;Generic_rD2]'*A)';


%% -------------------------------------------------------------------------
% Adjusted segment parameters Q and corresponding geometry
% -------------------------------------------------------------------------

% Pelvis
Segment(6).Informed.Q(4:6,:,:) = Segment(6).Informed.Q(4:6,:,:) + ...
    Mprod_array3(Segment(6).Informed.T(1:3,1:3,:), ...
    repmat([Scaled_rP6 - Informed_rP6],[1,1,n]));
Segment(6).Informed.Q(7:9,:,:) = Segment(6).Informed.Q(7:9,:,:) + ...
    Mprod_array3(Segment(6).Informed.T(1:3,1:3,:), ...
    repmat([Scaled_rD6 - Informed_rD6],[1,1,n]));

% Thigh
Segment(5).Informed.Q(4:6,:,:) = Segment(5).Informed.Q(4:6,:,:) + ...
    Mprod_array3(Segment(5).Informed.T(1:3,1:3,:), ...
    repmat([Scaled_rP5 - Informed_rP5],[1,1,n]));
Segment(5).Informed.Q(7:9,:,:) = Segment(5).Informed.Q(7:9,:,:) + ...
    Mprod_array3(Segment(5).Informed.T(1:3,1:3,:), ...
    repmat([Scaled_rD5 - Informed_rD5],[1,1,n]));

% Shank
Segment(3).Informed.Q(4:6,:,:) = Segment(5).Informed.Q(7:9,:,:); % rP3 = rD5
Segment(3).Informed.Q(7:9,:,:) = Segment(3).Informed.Q(7:9,:,:) + ...
    Mprod_array3(Segment(3).Informed.T(1:3,1:3,:), ...
    repmat([Scaled_rD3 - Informed_rD3],[1,1,n])); % rP3 = rD5

% Foot
Segment(2).Informed.Q(4:6,:,:) = Segment(3).Informed.Q(7:9,:,:); % rP2 = rD3
Segment(2).Informed.Q(7:9,:,:) = Segment(2).Informed.Q(7:9,:,:) + ...
    Mprod_array3(Segment(2).Informed.T(1:3,1:3,:), ...
    repmat([Scaled_rD2 - Informed_rD2],[1,1,n])); % rP2 = rD3


% Corresponding geometry
Segment(4).Informed.Scale = Model.Informed.Scale; % Patella not already defined
%
for i = [2,3,5,6] % From i = 2 (Foot) to i = 6 (Pelvis)
    
    % Scale
    Segment(i).Informed.Scale = Model.Informed.Scale;
    % Segment length
    Segment(i).Informed.L = Segment(i).Generic.L*Segment(i).Informed.Scale;
    % Segment angles
    Segment(i).Informed.alpha = Segment(i).Generic.alpha;
    Segment(i).Informed.beta = Segment(i).Generic.beta;
    Segment(i).Informed.gamma = Segment(i).Generic.gamma;

    % Matrix B from SCS to NSCS
    Segment(i).Informed.B = [1, ...
        Segment(i).Informed.L*cosd(Segment(i).Informed.gamma), ...
        cosd(Segment(i).Informed.beta); ...
        0, ...
        Segment(i).Informed.L*sind(Segment(i).Informed.gamma), ...
        (cosd(Segment(i).Informed.alpha) - ...
        cosd(Segment(i).Informed.beta)*cosd(Segment(i).Informed.gamma))/ ...
        sind(Segment(i).Informed.gamma); ...
        0, ...
        0, ...
        sqrt(1 - cosd(Segment(i).Informed.beta)^2 - ...
        ((cosd(Segment(i).Informed.alpha) - ...
        cosd(Segment(i).Informed.beta)*cosd(Segment(i).Informed.gamma))/ ...
        sind(Segment(i).Informed.gamma))^2)];
    
        % Transformation from ICS to SCS
        Segment(i).Informed.T = Q2Tuv_array3(Segment(i).Informed.Q);

end


function SSE = SSE(x,Xi,Xj)
A(2:4,1:3) = eye(3)*x(1,1); % Scale
A(1,1) = x(2,1); % Anterior-posterior translation
A(1,2) = x(3,1); % Superior-inferior translation
% Final control points minus transformed initial control points
Xja =  Xj - [ones(size(Xi,1),1),Xi]*A;
% Sum of squared errors (residual)
SSE = sqrt(sum(diag(Xja'*Xja)').^2); % RMS of 3 values on x,y,z

