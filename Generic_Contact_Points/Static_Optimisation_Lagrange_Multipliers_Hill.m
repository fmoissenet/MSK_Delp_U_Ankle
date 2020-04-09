% FUNCTION
% Static_Optimisation_Lagrange_Multipliers_Hill.m
%__________________________________________________________________________
%
% PURPOSE
% Computation musculo-tendon, contact, ligament and bone forces
%
% SYNOPSIS
% Model = Static_Optimisation_Lagrange_Multipliers_Hill(Segment,Joint,Model)
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
% Find the minimum of the sum of squared (cube) activations suject to equality
% (i.e., dynamic equlibrium) and inequality constraints (i.e., positive
% forces)
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
% Matlab R2020a
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
% fmincon algorithm = sqp (other than active-set)
%
% Modified by Raphael Dumas
% November 2019
% Activation-based criterion
%
% Modified by Raphael Dumas
% March 2020
% Generic vs. Informed structures (Segment, Joint, Model)
% n,f,fc as Model.Informed fields
%__________________________________________________________________________

function Model = Static_Optimisation_Lagrange_Multipliers_Hill(Segment,Joint,Model)

% Number of frames
n = Model.Informed.n;
% Acquisition frequency
f = Model.Informed.f;
% Cut off frequency for filtering
fc = Model.Informed.fc;

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
Model.Informed.K1(1,:,:) = Model.Informed.K(1,1:48,:) ... % Ankle contact about X axis
    *(-1); % Coefficient for contact force;
Model.Informed.K1(2,:,:) = Model.Informed.K(2,1:48,:) ... % Ankle contact about Y axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(3,:,:) = Model.Informed.K(3,1:48,:) ... % Ankle contact about Z axis
    *(-1); % Coefficient for contact force;
Model.Informed.K1(4,:,:) = Model.Informed.K(5,1:48,:) ... % Knee medial contact about X axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(5,:,:) = Model.Informed.K(6,1:48,:) ... % Knee medial contact about Y axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(6,:,:) = Model.Informed.K(7,1:48,:) ... % Knee medial contact about Z axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(7,:,:) = Model.Informed.K(8,1:48,:) ... % Knee lateral contact about X axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(8,:,:) = Model.Informed.K(9,1:48,:) ... % Knee lateral contact about Y axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(9,:,:) = Model.Informed.K(10,1:48,:) ... % Patellar contact about X axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(10,:,:) = Model.Informed.K(11,1:48,:) ... % Patellar contact about Y axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(11,:,:) = Model.Informed.K(12,1:48,:) ... % Patellar contact about Z axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(12,:,:) = Model.Informed.K(15,1:48,:) ... % PT
    *(1/(2*Joint(4).Informed.d(1,1))); % Coefficient for ligament force
Model.Informed.K1(13,:,:) = Model.Informed.K(16,1:48,:) ... % Hip contact about X axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(14,:,:) = Model.Informed.K(17,1:48,:) ... % Hip contact about Y axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(15,:,:) = Model.Informed.K(18,1:48,:) ...% Hip contact about Z axis
    *(-1); % Coefficient for contact force
Model.Informed.K1(16,:,:) = Model.Informed.K(22,1:48,:) ... % Foot axial
    *(-1/(2*Segment(2).Informed.L)); % Coefficient for bone force
Model.Informed.K1(17,:,:) = Model.Informed.K(28,1:48,:) ... % Tibia axial
    *(-1/(2*Segment(3).Informed.L)); % Coefficient for bone force
Model.Informed.K1(18,:,:) = Model.Informed.K(34,1:48,:) ... % Patella axial
    *(-1/(2*Segment(4).Informed.L)); % Coefficient for bone force
Model.Informed.K1(19,:,:) = Model.Informed.K(40,1:48,:) ... % Femur axial
    *(-1/(2*Segment(5).Informed.L)); % Coefficient for bone force

% Discarded forces
Model.Informed.K2(1,:,:) = Model.Informed.K(4,1:48,:); % Ankle angle
Model.Informed.K2(2,:,:) = Model.Informed.K(13,1:48,:); % Patella angle 1
Model.Informed.K2(3,:,:) = Model.Informed.K(14,1:48,:); % Patella angle 2
Model.Informed.K2(4,:,:) = Model.Informed.K(19,1:48,:);
Model.Informed.K2(5,:,:) = Model.Informed.K(20,1:48,:);
Model.Informed.K2(6,:,:) = Model.Informed.K(21,1:48,:);
Model.Informed.K2(7,:,:) = Model.Informed.K(23,1:48,:);
Model.Informed.K2(8,:,:) = Model.Informed.K(24,1:48,:);
Model.Informed.K2(9,:,:) = Model.Informed.K(25,1:48,:);
Model.Informed.K2(10,:,:) = Model.Informed.K(26,1:48,:);
Model.Informed.K2(11,:,:) = Model.Informed.K(27,1:48,:);
Model.Informed.K2(12,:,:) = Model.Informed.K(29,1:48,:);
Model.Informed.K2(13,:,:) = Model.Informed.K(30,1:48,:);
Model.Informed.K2(14,:,:) = Model.Informed.K(31,1:48,:);
Model.Informed.K2(15,:,:) = Model.Informed.K(32,1:48,:);
Model.Informed.K2(16,:,:) = Model.Informed.K(33,1:48,:);
Model.Informed.K2(17,:,:) = Model.Informed.K(35,1:48,:);
Model.Informed.K2(18,:,:) = Model.Informed.K(36,1:48,:);
Model.Informed.K2(19,:,:) = Model.Informed.K(37,1:48,:);
Model.Informed.K2(20,:,:) = Model.Informed.K(38,1:48,:);
Model.Informed.K2(21,:,:) = Model.Informed.K(39,1:48,:);
Model.Informed.K2(22,:,:) = Model.Informed.K(41,1:48,:);
Model.Informed.K2(23,:,:) = Model.Informed.K(42,1:48,:);

% -------------------------------------------------------------------------
% HILL MODEL
% -------------------------------------------------------------------------
% * Fla: Winter & Challis 2010
% * Flp: Schutte et al. 1992 (in Buchanan et al. 2004)
% * Fv : Selk Ghafari et al. 2009
% -------------------------------------------------------------------------

% Activation
for j = 1:43
    
    % Muscle length
    Model.Informed.Lm(j,1,1:n) = Mprod_array3((Model.Informed.Lmt(j,1,:) - ...
        Model.Informed.Lts(j)),repmat(cos(Model.Informed.pennation(j))^(-1),[1,1,n]));
    
    % Operaring range : 0.44*L0 (Fla = 0) to 1.2*L0 (Fpl too high)
    % Compliant tendon shall prevent too short/long muscle length
    Lmmin = 0.44*Model.Informed.L0(j);
    Lmmax = 1.2*Model.Informed.L0(j);
    Lmjmin = min(Model.Informed.Lm(j,1,1:n));
    if Lmjmin > 0.44*Model.Informed.L0(j)
        Lmmin = Lmjmin; % Unchanged for the muscle
    end
    Lmjmax = max(Model.Informed.Lm(j,1,1:n));
    if Lmjmax < 1.2*Model.Informed.L0(j)
        Lmmax = Lmjmax; % Unchanged for the muscle
    end
    % Rescaled centered on L0
    Model.Informed.Lm(j,1,1:n) = (Model.Informed.Lm(j,1,1:n) - ...
        repmat(Model.Informed.L0(j),[1,1,n])) *...
        (Lmmax-Lmmin)/(Lmjmax-Lmjmin) + repmat(Model.Informed.L0(j),[1,1,n]);
    % Translated if still Lm < 0.44*L0
    if min(Model.Informed.Lm(j,1,:)) <= 0.44*Model.Informed.L0(j)
        Model.Informed.Lm(j,1,:) = Model.Informed.Lm(j,1,:) + ...
            repmat((0.44*Model.Informed.L0(j) - ...
            min(Model.Informed.Lm(j,1,:)))*1.05,... 5% more
            [1,1,n]);
        display(['Adjusted length for muscle ',num2str(j)])
    end
    
    % Normalised muscle velocity
    Model.Informed.V(j,1,1:n) = Mprod_array3(Vfilt_array3(...
        Derive_array3(Model.Informed.Lm(j,:,:),1/f),f,fc), ...
        repmat(cos(Model.Informed.pennation(j)),[1,1,n]))./...
        repmat(Model.Informed.V0max(j),[1,1,n]);
    
    % Force-length
    Fla(j,1,1:n) = ones(1,1,n) - ...
        ((Model.Informed.Lm(j,1,:) - repmat(Model.Informed.L0(j),[1,1,n]))./...
        (0.56*repmat(Model.Informed.L0(j),[1,1,n]))).^2;
    
    % Force-velocity
    Fv(j,1,1:n) = ones(1,1,n) - tanh(0.3*Model.Informed.V(j,:,:));
    
    % Passive
    Flp(j,1,1:n) = exp(repmat(10,[1,1,n]).* ...
        (Model.Informed.Lm(j,1,:)./repmat(Model.Informed.L0(j),[1,1,n]) ... % Normalised
        - ones(1,1,n))).*repmat(1/exp(5),[1,1,n]);
    
%     % Less stiff musculo-tendon unit in older adults
%     Flp(j,1,1:n) = Flp(j,1,1:n)*0.70; % Less stiff musculo-tendon unit in older adults
    
    % From musculo-tendon force to activation
    % Fm = Fmax*cos(pennation) * (Fa(l)*F(v) * a + Fp(l))
    %  a = Fm/[Fa(l)*F(v)*Fmax*cos(pennation)] - Fp(l)/[Fa(l)*F(v)]
    % a = A*Fm + b
    A(j,j,1:n) = (repmat(61*Model.Informed.PCSA(j,1)*...
        cos(Model.Informed.pennation(j,1)),[1,1,n]).*...
        Fla(j,1,1:n).*Fv(j,1,1:n)).^-1;
    B(j,1,1:n) = - Flp(j,1,1:n)./(Fla(j,1,1:n).*Fv(j,1,1:n));
    
    % Mininal musculo-tendon force at activation = 0
    Fmin(j,1,1:n) = repmat(61*Model.Informed.PCSA(j,1)*...
        cos(Model.Informed.pennation(j,1)),[1,1,n]).*...
        Flp(j,1,1:n);
    
    % Maximal musculo-tendon force at activation = 1
    Fmax(j,1,1:n) = repmat(61*Model.Informed.PCSA(j,1)*...
        cos(Model.Informed.pennation(j,1)),[1,1,n]).*...
        (Fla(j,1,1:n).*Fv(j,1,1:n).*ones(1,1,n) + Flp(j,1,1:n));
    
end

Model.Informed.A = A; % To further transform f into a
Model.Informed.B = B; % To further transform f into a

% -------------------------------------------------------------------------
% RUN OPTIMISATION FOR EACH FRAME
% -------------------------------------------------------------------------
for i = 1:n
    
    waitbar(i/n,h);
    
    % ---------------------------------------------------------------------
    % Partial parameter reduction and constraint equations
    % ---------------------------------------------------------------------
    [eigvector,~] = eig(Model.Informed.K2(:,:,i)'*Model.Informed.K2(:,:,i));
    Model.Informed.ZK2(:,:,i) = eigvector(:,1:25); % 25 first eigenvalues are 0
    Aeq = Model.Informed.ZK2(:,:,i)'*[Model.Informed.Lever(:,:,i), - ...
        Model.Informed.K1(:,:,i)'];
    Beq = Model.Informed.ZK2(:,:,i)'*...
        (Model.Informed.G(:,:,i)*Model.Informed.d2Qdt2(:,:,i) - ...
        Model.Informed.P(:,:,i) - Model.Informed.R(:,:,i));
    
    % Initial guess
    if i == 1
        Xini = zeros(43 + 19,1);
    else
        Xini = Model.Informed.X(:,1,i-1);
    end
    
    % Lower abounds
    Xmin = zeros(43 + 19,1);
    %
    Xmin(43+1,1) = -Inf; % Ankle contact about X axis
    % Ankle contact about Y axis
    Xmin(43+3,1) = -Inf; % Ankle contact about Z axis
    Xmin(43+4,1) = -Inf; % Knee medial contact about X axis
    % Knee medial contact about Y axis
    Xmin(43+6,1) = -Inf; % Knee medial contact about Z axis
    Xmin(43+7,1) = -Inf; % Knee lateral contact about X axis
    % Knee lateral contact about Y axis
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
    
    % Muscle bounds from Hill model
    Xmin(1:43) = Fmin(:,:,i);
    Xmax(1:43) = Fmax(:,:,i);
    
    
    % ---------------------------------------------------------------------
    % Run optimization (fmincon)
    % ---------------------------------------------------------------------
    options = optimset('Display','notify-detailed','LargeScale','off',...
        'MaxFunEvals',8000,...  Default is 100*numberOfVariables
        'MaxIter',800,... Default is 400
        'TolFun',1e-4,... Default is 1e-6
        'TolX',1e-6,... Default is 1e-6
        'TolCon',1e-6,... Default is 1e-6
        'GradObj','on','algorithm','sqp');
    
    % Optimisation
    Aw(1:43,1:43) = A(:,:,i); % From musculo-tendon force to activation
    Aw(43+1:43+19,43+1:43+19) = diag(zeros(19,1));
%         Aw(43+2,43+2) = 1/3000*43/4; % Ankle contact about Y axis
%         Aw(43+5,43+5) = 1/1500*43/4; % Knee medial contact about Y axis
%         Aw(43+8,43+8) = 1/1500*43/4; % Knee lateral contact about Y axis
%         Aw(43+14,43+14) = 1/3000*43/4; % Hip contact about Y axis
    Bw = [B(:,:,i); ... % From musculo-tendon force to activation
        zeros(19,1)];
    C = @(X)Criterion_Lagrange_Multipliers(X,Aw,Bw); % Anonymous function (to pass extra parameters)
    [X,~,exitflag] = fmincon(C,Xini,[],[],Aeq,Beq,Xmin,Xmax,[],options);
    Model.Informed.X(:,1,i) = X;
    flag(i) = exitflag;
    
end

close(h);
% figure
% plot(flag);

% -------------------------------------------------------------------------
% SUBFUNCTION
% -------------------------------------------------------------------------

function [J,G] = Criterion_Lagrange_Multipliers(X,A,B)

% Objective function
J = 1/2*(A*X+B)'*(A*X+B); % Sum of squared activations
% Gradient
G = A*(A*X+B);

