clearvars;
clc;

% Model landmarks
RASIS = [0.0; 0.0; 0.128]; % Directly in ground
LASIS = [0.0; 0.0; -0.128]; % Directly in ground
RPSIS = [-0.148; 0.0; 0.047]; % Directly in ground
LPSIS = [-0.148; 0.0; -0.047]; % Directly in ground
LJC   = [-0.1055; 0.072; 0]; % Directly in ground
RHJC  = [-0.0707; -0.0661; 0.0835]; % Directly in ground
mASIS = (RASIS+LASIS)/2; % mid ASIS

% Apply flesh margin (1 cm along X for ASIS landmarks to introduce a flesh 
% width)
RASIS = RASIS+[0.01;0;0];
LASIS = LASIS+[0.01;0;0];

% Adjusted pelvis SCS
Z_pelvis = (RASIS-LASIS)/norm(RASIS-LASIS);
Y_pelvis = cross(((RPSIS+LPSIS)/2-RASIS),((RPSIS+LPSIS)/2-LASIS))/ ...
           norm(cross(((RPSIS+LPSIS)/2-RASIS),((RPSIS+LPSIS)/2-LASIS)));
X_pelvis = cross(Y_pelvis,Z_pelvis)/norm(cross(Y_pelvis,Z_pelvis));
T_pelvis = [X_pelvis Y_pelvis Z_pelvis mASIS; 0 0 0 1];
W_pelvis = norm(RASIS-LASIS); % Pelvis width of the generic musculoskeletal model

% Set all landmarks in the pelvis adjusted SCS
temp  = T_pelvis*[RASIS;1];
RASIS = temp(1:3);
clear temp;
temp  = T_pelvis*[LASIS;1];
LASIS = temp(1:3);
clear temp;
temp  = T_pelvis*[RPSIS;1];
RPSIS = temp(1:3);
clear temp;
temp  = T_pelvis*[LPSIS;1];
LPSIS = temp(1:3);
clear temp;
temp  = T_pelvis*[LJC;1];
LJC = temp(1:3);
clear temp;
temp  = T_pelvis*[RHJC;1];
RHJC = temp(1:3);
clear temp;

% Compute regressions
RHJC_reg = RHJC/W_pelvis;
LHJC_reg = [RHJC_reg(1);RHJC_reg(2);-RHJC_reg(3)];
LJC_reg  = LJC/W_pelvis;

figure; 
hold on;
axis equal;
plot3(RASIS(1),RASIS(2),RASIS(3),'Marker','o','Color','red');
plot3(LASIS(1),LASIS(2),LASIS(3),'Marker','o','Color','red');
plot3(RPSIS(1),RPSIS(2),RPSIS(3),'Marker','o','Color','red');
plot3(LPSIS(1),LPSIS(2),LPSIS(3),'Marker','o','Color','red');
plot3(RHJC(1),RHJC(2),RHJC(3),'Marker','o','Color','blue');
plot3(LJC(1),LJC(2),LJC(3),'Marker','o','Color','green');
quiver3(mASIS(1),mASIS(2),mASIS(3),X_pelvis(1),X_pelvis(2),X_pelvis(3),0.1,'Color','red');
quiver3(mASIS(1),mASIS(2),mASIS(3),Y_pelvis(1),Y_pelvis(2),Y_pelvis(3),0.1,'Color','blue');
quiver3(mASIS(1),mASIS(2),mASIS(3),Z_pelvis(1),Z_pelvis(2),Z_pelvis(3),0.1,'Color','green');

save('HJC_LCJ_regression_Delp.mat','RHJC_reg','LHJC_reg','LJC_reg');
