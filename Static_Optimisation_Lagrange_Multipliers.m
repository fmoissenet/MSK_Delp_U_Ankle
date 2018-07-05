% FUNCTION
% Static_Optimisation_Lagrange_Multipliers.m
%__________________________________________________________________________
%
% PURPOSE
% Computation musculo-tendon, contact, ligament and bone forces
%
% SYNOPSIS
% Model = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Joint (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% OUTPUT
% Model (cf. data structure in user guide)
%
% DESCRIPTION
% Find the minimum of the sum of squared forces suject to equality 
% (i.e., dynamic equlibrium) and inequality constraints (i.e., positive
% forces
%
% REFERENCE
% F Moissenet, L Cheze, R Dumas. A 3D lower limb musculoskeletal model for
% simultaneous estimation of musculo-tendon, joint contact, ligament and
% bone forces during gait. Journal of Biomechanics 2014;47(1):50-8.
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
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
% February 2014
% Full selection of lambda1
%
% Modified by Raphael Dumas
% April 2015
% Coefficient and bounds adapted to contact point constraints
%
% Modified by Raphael Dumas
% January 2018
% Universal joint at the ankle
% Anonymous function (to pass extra parameters)
%
% Modified by Raphael Dumas
% April 2018
% fimoncon algorithm = sqp (other than active-set)
%__________________________________________________________________________
%
% Copyright (C) 2018  Raphael Dumas, Florent Moissenet
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%__________________________________________________________________________

function Model = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model,Force,weights,method)

% Number of frames
n = size(Model.K,3);

% Waitbar
h = waitbar(1,'Optimization is running...');

% -------------------------------------------------------------------------
% Prepare reduced Jacobian matrices
% K1 is kept, K2 is removed
%
% Coefficients are introduced in order to associate Lagrange multiplieurs
% with external forces acting on the proximal segment in order to have 
% positive force about the SCS axes (or axes of the structure)
% Contact forces: reaction forces
% Ligament forces: tensile forces
% Bone forces: compression forces
% -------------------------------------------------------------------------

% Included contact, ligament and bone forces
Model.K1(1,:,:) = Model.K(1,1:48,:) ... % Ankle contact about X axis
    *(-1); % Coefficient for contact force;
Model.K1(2,:,:) = Model.K(2,1:48,:) ... % Ankle contact about Y axis 
    *(-1); % Coefficient for contact force
Model.K1(3,:,:) = Model.K(3,1:48,:) ... % Ankle contact about Z axis
    *(-1); % Coefficient for contact force;
Model.K1(4,:,:) = Model.K(5,1:48,:) ... % Knee medial contact 
    *(-1); % Coefficient for contact force
Model.K1(5,:,:) = Model.K(6,1:48,:) ... % Knee lateral contact 
    *(-1); % Coefficient for contact force
Model.K1(6,:,:) = Model.K(7,1:48,:) ... % ACL 
    *(1/(2*Joint(3).d(1,3))); % Coefficient for ligament force
Model.K1(7,:,:) = Model.K(8,1:48,:) ... % PCL 
    *(1/(2*Joint(3).d(1,4))); % Coefficient for ligament force
Model.K1(8,:,:) = Model.K(9,1:48,:) ... % MCL
    *(1/(2*Joint(3).d(1,4))); % Coefficient for ligament force;
Model.K1(9,:,:) = Model.K(10,1:48,:) ... % Patellar contact about X axis
    *(-1); % Coefficient for contact force
Model.K1(10,:,:) = Model.K(11,1:48,:) ... % Patellar contact about Y axis
    *(-1); % Coefficient for contact force
Model.K1(11,:,:) = Model.K(12,1:48,:) ... % Patellar contact about Z axis
    *(-1); % Coefficient for contact force
Model.K1(12,:,:) = Model.K(15,1:48,:) ... % PT 
    *(1/(2*Joint(4).d(1,1))); % Coefficient for ligament force
Model.K1(13,:,:) = Model.K(16,1:48,:) ... % Hip contact about X axis 
    *(-1); % Coefficient for contact force
Model.K1(14,:,:) = Model.K(17,1:48,:) ... % Hip contact about Y axis 
    *(-1); % Coefficient for contact force
Model.K1(15,:,:) = Model.K(18,1:48,:) ...% Hip contact about Z axis 
    *(-1); % Coefficient for contact force
Model.K1(16,:,:) = Model.K(22,1:48,:) ... % Foot axial 
    *(-1/(2*Segment(2).L)); % Coefficient for bone force
Model.K1(17,:,:) = Model.K(28,1:48,:) ... % Tibia axial 
    *(-1/(2*Segment(3).L)); % Coefficient for bone force
Model.K1(18,:,:) = Model.K(34,1:48,:) ... % Patella axial 
    *(-1/(2*Segment(4).L)); % Coefficient for bone force
Model.K1(19,:,:) = Model.K(40,1:48,:) ... % Femur axial 
    *(-1/(2*Segment(5).L)); % Coefficient for bone force

% Discarded forces
Model.K2(1,:,:) = Model.K(4,1:48,:); % Ankle angle
Model.K2(2,:,:) = Model.K(13,1:48,:); 
Model.K2(3,:,:) = Model.K(14,1:48,:);
Model.K2(4,:,:) = Model.K(19,1:48,:);
Model.K2(5,:,:) = Model.K(20,1:48,:);
Model.K2(6,:,:) = Model.K(21,1:48,:);
Model.K2(7,:,:) = Model.K(23,1:48,:);
Model.K2(8,:,:) = Model.K(24,1:48,:);
Model.K2(9,:,:) = Model.K(25,1:48,:);
Model.K2(10,:,:) = Model.K(26,1:48,:);
Model.K2(11,:,:) = Model.K(27,1:48,:);
Model.K2(12,:,:) = Model.K(29,1:48,:);
Model.K2(13,:,:) = Model.K(30,1:48,:);
Model.K2(14,:,:) = Model.K(31,1:48,:);
Model.K2(15,:,:) = Model.K(32,1:48,:);
Model.K2(16,:,:) = Model.K(33,1:48,:);
Model.K2(17,:,:) = Model.K(35,1:48,:);
Model.K2(18,:,:) = Model.K(36,1:48,:);
Model.K2(19,:,:) = Model.K(37,1:48,:);
Model.K2(20,:,:) = Model.K(38,1:48,:);
Model.K2(21,:,:) = Model.K(39,1:48,:);
Model.K2(22,:,:) = Model.K(41,1:48,:);
Model.K2(23,:,:) = Model.K(42,1:48,:);

% -------------------------------------------------------------------------
% Weight matrix
% -------------------------------------------------------------------------
W = eye(43+19); % Number of muscles + number of Lagrange multipliers
for i = 1:43
    W(i,i) = 1;%(1/43);
end
for i = 1:19
    W(43+i,43+i) = weights(i);
end

% -------------------------------------------------------------------------
% RUN OPTIMISATION FOR EACH FRAME
% -------------------------------------------------------------------------
if strcmp(method,'mmm')
    for i = 1:n

%         waitbar(i/n,h);

        % ---------------------------------------------------------------------
        % Partial parameter reduction and constraint equations
        % ---------------------------------------------------------------------
        [eigvector,~] = eig(Model.K2(:,:,i)'*Model.K2(:,:,i));
        Model.ZK2(:,:,i) = eigvector(:,1:25); % 25 first eigenvalues are 0
        Aeq = Model.ZK2(:,:,i)'*[Model.Lever(:,:,i), - Model.K1(:,:,i)'];
        Beq = Model.ZK2(:,:,i)'*...
            (Model.G(:,:,i)*Model.d2Qdt2(:,:,i) - Model.P(:,:,i) - Model.R(:,:,i));

        % Initial guess
        if i == 1
            Xini = zeros(43 + 19,1);
        else
            Xini = Model.X(:,1,i-1);
        end

        % Lower abounds
        Xmin = zeros(43 + 19,1);
        Xmin(43+1,1) = -Inf; % Ankle contact about X axis
        % Ankle contact about Y axis
        Xmin(43+3,1) = -Inf; % Ankle contact about Z axis
        % Knee medial contact
        % Knee lateral contact
        % ACL
        % PCL
        Xmin(43+8,1) = -Inf; % MCL
        Xmin(43+9,1) = -Inf; % Patellar contact about X axis
        Xmin(43+10,1) = -Inf; % Patellar contact about Y axis
        Xmin(43+11,1) = -Inf; % Patellar contact about Z axis
        % PT
        Xmin(43+13,1) = -Inf; % Hip contact about X axis
        % Hip contact about Y axis
        Xmin(43+15,1) = -Inf; % Hip contact about Z axis
        Xmin(43+16,1) = -Inf; % Foot axial
        % Tibia axial
        Xmin(43+18,1) = -Inf; % Patella axial
        % Femur axial

        % Upper bounds
        Xmax = Inf(43 + 19,1);
    %     Xmax = [Model.Fmax(1:43,1);Inf(19,1)];
%         acoeff = 1.3;
%         Xmax = [Model.PCSA(1:43,1)*61;Inf(19,1)]*acoeff; % Muscle stress of 61 N/cm^2
%         Xmax(43+[4;5]) = [max(Force.KneeMedial);max(Force.KneeLateral)]*acoeff;
%         Xmax(43+[6;7]) = [3000;2000]*acoeff; % Cleather and Bull, 2011

        % ---------------------------------------------------------------------
        % Run optimization (fmincon)
        % ---------------------------------------------------------------------
        options = optimset('Display','off','MaxIter',500,'LargeScale','off','TolFun',0.1,'TolX',1e-10,'TolCon',1e-10,...
            'algorithm','sqp','GradObj','off');
        % Weighted sum method
        C = @(X)Criterion_Lagrange_Multipliers_WeightedSum(X,W,43); % Anonymous function (to pass extra parameters)
        [X,~,exitflag] = fmincon(C,Xini,[],[],Aeq,Beq,Xmin,Xmax,[],options);
        Model.X(:,1,i) = X;
        flag(i) = exitflag;
            
        m = 43;
        clear x;
        x = 1:m;
        Jtest(1,:,i) = ((X(x,1))'*(X(x,1)));
        x = m+[4;5]; % only tibiofemoral contact forces
        Jtest(2,:,i) = ((X(x,1))'*(X(x,1)));
        x = m+[6;7]; % only cruciate ligaments
        Jtest(3,:,i) = ((X(x,1))'*(X(x,1)));

    end
%     figure
% plot(flag);
%     
    Jmax(1) = max(Jtest(1,:,:)); 
    Jmax(2) = max(Jtest(2,:,:)); 
    Jmax(3) = max(Jtest(3,:,:)); 
%     
    figure;
hold on;
plot(squeeze(Jtest(1,:,:))','red');
plot(squeeze(Jtest(2,:,:))','green');
plot(squeeze(Jtest(3,:,:))','blue');    
end

for i = 1:n
    
    waitbar(i/n,h);
    
    % ---------------------------------------------------------------------
    % Partial parameter reduction and constraint equations
    % ---------------------------------------------------------------------
    [eigvector,~] = eig(Model.K2(:,:,i)'*Model.K2(:,:,i));
    Model.ZK2(:,:,i) = eigvector(:,1:25); % 25 first eigenvalues are 0
    Aeq = Model.ZK2(:,:,i)'*[Model.Lever(:,:,i), - Model.K1(:,:,i)'];
    Beq = Model.ZK2(:,:,i)'*...
        (Model.G(:,:,i)*Model.d2Qdt2(:,:,i) - Model.P(:,:,i) - Model.R(:,:,i));
    
    % Initial guess
    if i == 1
        Xini = zeros(43 + 19,1);
    else
        Xini = Model.X(:,1,i-1);
    end
    
    % Lower abounds
    Xmin = zeros(43 + 19,1);
    Xmin(43+1,1) = -Inf; % Ankle contact about X axis
    % Ankle contact about Y axis
    Xmin(43+3,1) = -Inf; % Ankle contact about Z axis
    % Knee medial contact
    % Knee lateral contact
    % ACL
    % PCL
    Xmin(43+8,1) = -Inf; % MCL
    Xmin(43+9,1) = -Inf; % Patellar contact about X axis
    Xmin(43+10,1) = -Inf; % Patellar contact about Y axis
    Xmin(43+11,1) = -Inf; % Patellar contact about Z axis
    % PT
    Xmin(43+13,1) = -Inf; % Hip contact about X axis
    % Hip contact about Y axis
    Xmin(43+15,1) = -Inf; % Hip contact about Z axis
    Xmin(43+16,1) = -Inf; % Foot axial
    % Tibia axial
    Xmin(43+18,1) = -Inf; % Patella axial
    % Femur axial
    
    % Upper bounds
    Xmax = Inf(43 + 19,1);
%     Xmax = [Model.Fmax(1:43,1);Inf(19,1)];
%     acoeff = 1.3;
%     Xmax = [Model.PCSA(1:43,1)*61;Inf(19,1)]*acoeff; % Muscle stress of 61 N/cm^2
%     Xmax(43+[4;5]) = [max(Force.KneeMedial);max(Force.KneeLateral)]*acoeff;
%     Xmax(43+[6;7]) = [3000;2000]*acoeff; % Cleather and Bull, 2011
    
    % ---------------------------------------------------------------------
    % Run optimization (fmincon)
    % ---------------------------------------------------------------------
    options = optimset('Display','off','MaxIter',500,'LargeScale','off','TolFun',0.1,'TolX',1e-10,'TolCon',1e-10,...
        'algorithm','sqp','GradObj','off');
    if strcmp(method,'wsm')
        % Weighted sum method
        C = @(X)Criterion_Lagrange_Multipliers_WeightedSum(X,W,43); % Anonymous function (to pass extra parameters)
        [X,J,exitflag] = fmincon(C,Xini,[],[],Aeq,Beq,Xmin,Xmax,[],options);
    elseif strcmp(method,'mmm')
        % Min max method
        C = @(X)Criterion_Lagrange_Multipliers_MinMax(X,Jmax,43); % Anonymous function (to pass extra parameters)
        [X,J,~,exitflag] = fminimax(C,Xini,[],[],Aeq,Beq,Xmin,Xmax,[],options);
    end

    Model.X(:,1,i) = X;
    Model.J(:,1,i) = J;
    flag(i) = exitflag;
    
end

close(h);
figure
plot(flag);

figure;
hold on;
plot(squeeze(Model.J(1,:,:))','red');
plot(squeeze(Model.J(2,:,:))','green');
plot(squeeze(Model.J(3,:,:))','blue');

% GC1, cycle 1 : max J2 at 34 and 97
% = Bosser en particulier sur ces deux points
% = produire le front de pareto avec  gamultiobj
% = tracer qqch de similaire que http://www.ozeninc.com/wp-content/uploads/2015/07/methods3.gif
% = tracer le plan (?) correspondant aux mesures (sur J2 a priori seulement)
% = positionner WSM1, WSM2 et MMM sur ce front de Pareto
% = est-ce que ces solutions sont bien sur le front?
% = est-ce qu'il y en a une qui se rapproche du plan de mesure?

% -------------------------------------------------------------------------
% SUBFUNCTIONS
% -------------------------------------------------------------------------
function J = Criterion_Lagrange_Multipliers_WeightedSum(X,W,m)

% J1: musculo-tendon forces
% J2: contact forces
% J3: ligament forces

x = 1:m;
J1 = ((X(x,1))'*W(x,x)*(X(x,1)));

x = m+[4;5]; % only tibiofemoral contact forces
J2 = ((X(x,1))'*W(x,x)*(X(x,1)));

x = m+[6;7]; % only cruciate ligaments
J3 = ((X(x,1))'*W(x,x)*(X(x,1)));

J = J1+J2+J3;

function J = Criterion_Lagrange_Multipliers_MinMax(X,Jmax,m)

% J1: musculo-tendon forces
% J2: contact forces
% J3: ligament forces

x = 1:m;
J(1) = (1/Jmax(1))*((X(x,1))'*(X(x,1)));

x = m+[4;5]; % only tibiofemoral contact forces
J(2) = (1/Jmax(2))*((X(x,1))'*(X(x,1)));

x = m+[6;7]; % only cruciate ligaments
J(3) = (1/Jmax(3))*((X(x,1))'*(X(x,1)));