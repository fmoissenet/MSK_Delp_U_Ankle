% Q2DoFs
function Dofs = Q2DoFs(Q)

% Segment
% Q(1:12,1);     % Foot
% Q(13:24,1);    % Shank
% Q(25:36,1);    % Patella
% Q(37:48,1);    % Thigh
% Q(49:60,1);    % Pelvis

% Ankle
R2 = Q2Rwu(Q(13:24,1))\Q2Ruw(Q(1:12,1));
Euler2 = R2mobileZYX_array3(R2);
Dofs(1,1) = Euler2(1,1); % Flexion-Extension (about Z=w proximal SCS axis)
Dofs(2,1) = Euler2(1,3); % Adduction-Abduction (about X=u distal SCS axis)

% Knee
R3 = Q2Rwu(Q(37:48,1))\Q2Ruv(Q(13:24,1));
Euler3 = R2mobileZXY_array3(R3); 
Dofs(3,1) = Euler3(1,1); % Flexion-Extension (about Z=w proximal SCS axis)

% Hip
R5 = Q2Rwu(Q(49:60,1))\Q2Ruv(Q(37:48,1));
Euler5 = R2mobileZXY_array3(R5);
Dofs(4,1) = Euler5(1,1); % Flexion-Extension (about Z=w proximal SCS axis)
Dofs(5,1) = Euler5(1,2); % Adduction-Abduction (about floatting axis)
Dofs(6,1) = Euler5(1,3); % Internal-External Rotation (about Y distal SCS axis)

end

% Subfunctions
function R = Q2Ruw(Q)
X = Q(1:3,1,:); % X = u
Y = cross(Q(10:12,1,:),X); % w x (u)
Y = Y/norm(Y);
Z = cross(X,Y);
Z = Z/norm(Z);
R = [X, Y, Z]; % Rotation matrix
end

function R = Q2Ruv(Q)
X = Q(1:3,1,:); % X = u
Z = cross(X,Q(4:6,1,:) - Q(7:9,1,:)); % X x (rP - rD)
Z = Z/norm(Z);
Y = cross(Z,X);
Y = Y/norm(Y);
R = [X, Y, Z]; % Rotation matrix
end

function R = Q2Rwu(Q)
Z = Q(10:12,1,:); % Z = w
Y = cross(Z,Q(1:3,1,:)); % Z x u
Y = Y/norm(Y);
X = cross(Y,Z);
X = X/norm(X);
R = [X, Y, Z]; % Rotation matrix
end

