% MAIN PROGRAM
% Main_Model_Visualisation.m
%__________________________________________________________________________
%
% PURPOSE
% Plotting of segment, joint and muscle features
%
% SYNOPSIS
% N/A (i.e., main program)
%
% DESCRIPTION
% Plotting of segment axes (u, w), endpoints (rP, rD) and markers (rM)
% (cf. data structure in user guide), of contacts and ligaments, and of
% muscle lines of actions at a set of frames
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Mprod_array3.m
% Q2Tuv_array3.m
% 
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphael Dumas
% September 2012
%
% Modified by Raphael Dumas
% April 2018
% Q2Tuv other than Q2Tuw to compute the positions of the centres of mass
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


% Figure
figure;
hold on;
axis equal;

% ICS
quiver3(0,0,0,1,0,0,0.5,'k');
quiver3(0,0,0,0,1,0,0.5,'k');
quiver3(0,0,0,0,0,1,0.5,'k');

% Number of frames
n = size(Segment(2).Q,3); % At least Q parameter for foot (or hand) segment

% Frames of interest
ni = [1,round(n*60/100),n]; 

%% Segments

% Forceplate
% Proximal enpoint P: P1 = centre of pressure
P1x = permute(Segment(1).Q(4,1,ni),[3,2,1]);
P1y = permute(Segment(1).Q(5,1,ni),[3,2,1]);
P1z = permute(Segment(1).Q(6,1,ni),[3,2,1]);
plot3(P1x,P1y,P1z,'ok');


% Foot
% Proximal enpoint P
P2x = permute(Segment(2).Q(4,1,ni),[3,2,1]);
P2y = permute(Segment(2).Q(5,1,ni),[3,2,1]);
P2z = permute(Segment(2).Q(6,1,ni),[3,2,1]);
plot3(P2x,P2y,P2z,'ob' );
% Distal endpoints D
D2x = permute(Segment(2).Q(7,1,ni),[3,2,1]);
D2y = permute(Segment(2).Q(8,1,ni),[3,2,1]);
D2z = permute(Segment(2).Q(9,1,ni),[3,2,1]);
plot3(D2x,D2y,D2z,'.b' );
plot3([P2x,D2x]',[P2y,D2y]',[P2z,D2z]','b');
% U axis
U2x = permute(Segment(2).Q(1,1,ni), [3,2,1]);
U2y = permute(Segment(2).Q(2,1,ni), [3,2,1]);
U2z = permute(Segment(2).Q(3,1,ni), [3,2,1]);
quiver3(P2x,P2y,P2z,U2x,U2y,U2z,0.25,'b');
% W axis
W2x = permute(Segment(2).Q(10,1,ni),[3,2,1]);
W2y = permute(Segment(2).Q(11,1,ni),[3,2,1]);
W2z = permute(Segment(2).Q(12,1,ni),[3,2,1]);
quiver3(D2x,D2y,D2z,W2x,W2y,W2z,0.25,'b');
% Markers
for j = 1:size(Segment(2).rM,2)
    plot3(permute(Segment(2).rM(1,j,ni),[3,2,1]),...
        permute(Segment(2).rM(2,j,ni),[3,2,1]),...
        permute(Segment(2).rM(3,j,ni),[3,2,1]),'+b');
end
% Centre of mass
G2 = Mprod_array3(Q2Tuv_array3(Segment(2).Q),...
    repmat([Segment(2).rCs;1],[1 1 n]));
G2x = permute(G2(1,1,ni),[3,2,1]);
G2y = permute(G2(2,1,ni),[3,2,1]);
G2z = permute(G2(3,1,ni),[3,2,1]);
plot3(G2x,G2y,G2z,'*b');


% Shank
% Proximal enpoint P
P3x = permute(Segment(3).Q(4,1,ni),[3,2,1]);
P3y = permute(Segment(3).Q(5,1,ni),[3,2,1]);
P3z = permute(Segment(3).Q(6,1,ni),[3,2,1]);
plot3(P3x,P3y,P3z,'or');
% Distal endpoints D
D3x = permute(Segment(3).Q(7,1,ni), [3,2,1]);
D3y = permute(Segment(3).Q(8,1,ni), [3,2,1]);
D3z = permute(Segment(3).Q(9,1,ni), [3,2,1]);
plot3(D3x,D3y,D3z,'.r');
plot3([P3x,D3x]',[P3y,D3y]',[P3z,D3z]','r');
% U axis
U3x = permute(Segment(3).Q(1,1,ni), [3,2,1]);
U3y = permute(Segment(3).Q(2,1,ni), [3,2,1]);
U3z = permute(Segment(3).Q(3,1,ni), [3,2,1]);
quiver3(P3x,P3y,P3z,U3x,U3y,U3z,0.25,'r');
% W axis
W3x = permute(Segment(3).Q(10,1,ni), [3,2,1]);
W3y = permute(Segment(3).Q(11,1,ni), [3,2,1]);
W3z = permute(Segment(3).Q(12,1,ni), [3,2,1]);
quiver3(D3x,D3y,D3z,W3x,W3y,W3z,0.25,'r');
% Markers
for j = 1:size(Segment(3).rM,2)
    plot3(permute(Segment(3).rM(1,j,ni),[3,2,1]),...
        permute(Segment(3).rM(2,j,ni),[3,2,1]),...
        permute(Segment(3).rM(3,j,ni),[3,2,1]),'+r');
end
% Centre of mass
G3 = Mprod_array3(Q2Tuv_array3(Segment(3).Q),...
    repmat([Segment(3).rCs;1],[1 1 n]));
G3x = permute(G3(1,1,ni),[3,2,1]);
G3y = permute(G3(2,1,ni),[3,2,1]);
G3z = permute(G3(3,1,ni),[3,2,1]);
plot3(G3x,G3y,G3z,'*r');


% Patella
% Proximal enpoint P
P4x = permute(Segment(4).Q(4,1,ni),[3,2,1]);
P4y = permute(Segment(4).Q(5,1,ni),[3,2,1]);
P4z = permute(Segment(4).Q(6,1,ni),[3,2,1]);
plot3(P4x,P4y,P4z,'og');
% Distal endpoints D
D4x = permute(Segment(4).Q(7,1,ni), [3,2,1]);
D4y = permute(Segment(4).Q(8,1,ni), [3,2,1]);
D4z = permute(Segment(4).Q(9,1,ni), [3,2,1]);
plot3(D4x,D4y,D4z,'.g');
plot3([P4x,D4x]',[P4y,D4y]',[P4z,D4z]','g');
% U axis
U4x = permute(Segment(4).Q(1,1,ni), [3,2,1]);
U4y = permute(Segment(4).Q(2,1,ni), [3,2,1]);
U4z = permute(Segment(4).Q(3,1,ni), [3,2,1]);
quiver3(P4x,P4y,P4z,U4x,U4y,U4z,0.25,'g');
% W axis
W4x = permute(Segment(4).Q(10,1,ni), [3,2,1]);
W4y = permute(Segment(4).Q(11,1,ni), [3,2,1]);
W4z = permute(Segment(4).Q(12,1,ni), [3,2,1]);
quiver3(D4x,D4y,D4z,W4x,W4y,W4z,0.25,'g');
% No markers
G4 = Mprod_array3(Q2Tuv_array3(Segment(4).Q),...
    repmat([Segment(4).rCs;1],[1 1 n]));
G4x = permute(G4(1,1,ni),[3,2,1]);
G4y = permute(G4(2,1,ni),[3,2,1]);
G4z = permute(G4(3,1,ni),[3,2,1]);
plot3(G4x,G4y,G4z,'*g'); % Should be supperimposed with P4

% Thigh
% Proximal enpoint P
P5x = permute(Segment(5).Q(4,1,ni), [3,2,1]);
P5y = permute(Segment(5).Q(5,1,ni), [3,2,1]);
P5z = permute(Segment(5).Q(6,1,ni), [3,2,1]);
plot3(P5x,P5y,P5z,'oc');
% Distal endpoints D
D5x = permute(Segment(5).Q(7,1,ni), [3,2,1]);
D5y = permute(Segment(5).Q(8,1,ni), [3,2,1]);
D5z = permute(Segment(5).Q(9,1,ni), [3,2,1]);
plot3(D5x,D5y,D5z,'.c');
plot3([P5x,D5x]',[P5y,D5y]',[P5z,D5z]','c');
% U axis
U5x = permute(Segment(5).Q(1,1,ni), [3,2,1]);
U5y = permute(Segment(5).Q(2,1,ni), [3,2,1]);
U5z = permute(Segment(5).Q(3,1,ni), [3,2,1]);
quiver3(P5x,P5y,P5z,U5x,U5y,U5z,0.25,'c');
% W axis
W5x = permute(Segment(5).Q(10,1,ni), [3,2,1]);
W5y = permute(Segment(5).Q(11,1,ni), [3,2,1]);
W5z = permute(Segment(5).Q(12,1,ni), [3,2,1]);
quiver3(D5x,D5y,D5z,W5x,W5y,W5z,0.25,'c');
% Markers
for j = 1:size(Segment(5).rM,2)
    plot3(permute(Segment(5).rM(1,j,ni),[3,2,1]),...
        permute(Segment(5).rM(2,j,ni),[3,2,1]),...
        permute(Segment(5).rM(3,j,ni),[3,2,1]),'+c');
end
% Centre of mass
G5 = Mprod_array3(Q2Tuv_array3(Segment(5).Q),...
    repmat([Segment(5).rCs;1],[1 1 n]));
G5x = permute(G5(1,1,ni),[3,2,1]);
G5y = permute(G5(2,1,ni),[3,2,1]);
G5z = permute(G5(3,1,ni),[3,2,1]);
plot3(G5x,G5y,G5z,'*c');


% Pelvis
% Proximal enpoint P
P6x = permute(Segment(6).Q(4,1,ni), [3,2,1]);
P6y = permute(Segment(6).Q(5,1,ni), [3,2,1]);
P6z = permute(Segment(6).Q(6,1,ni), [3,2,1]);
plot3(P6x,P6y,P6z,'om');
% Distal endpoints D
D6x = permute(Segment(6).Q(7,1,ni), [3,2,1]);
D6y = permute(Segment(6).Q(8,1,ni), [3,2,1]);
D6z = permute(Segment(6).Q(9,1,ni), [3,2,1]);
plot3(D6x,D6y,D6z,'.m');
plot3([P6x,D6x]',[P6y,D6y]',[P6z,D6z]','m');
% U axis
U6x = permute(Segment(6).Q(1,1,ni), [3,2,1]);
U6y = permute(Segment(6).Q(2,1,ni), [3,2,1]);
U6z = permute(Segment(6).Q(3,1,ni), [3,2,1]);
quiver3(P6x,P6y,P6z,U6x,U6y,U6z,0.25,'m');
% W axis
W6x = permute(Segment(6).Q(10,1,ni), [3,2,1]);
W6y = permute(Segment(6).Q(11,1,ni), [3,2,1]);
W6z = permute(Segment(6).Q(12,1,ni), [3,2,1]);
quiver3(D6x,D6y,D6z,W6x,W6y,W6z,0.25,'m');
% Markers
for j = 1:size(Segment(6).rM,2)
    plot3(permute(Segment(6).rM(1,j,ni),[3,2,1]),...
        permute(Segment(6).rM(2,j,ni),[3,2,1]),...
        permute(Segment(6).rM(3,j,ni),[3,2,1]),'+m');
end
% No center of mass


%% Joints

% Ankle
% Subtalar axis position in foot segment
% V12: virtual marker 1 of segment 2
V12 = Mprod_array3(repmat(([Segment(2).nV(1,1)*eye(3),...
    (1 + Segment(2).nV(2,1))*eye(3), ...
    - Segment(2).nV(2,1)*eye(3), ...
    Segment(2).nV(3,1)*eye(3)]),[1 1 n]),Segment(2).Q);  
V12x = permute(V12(1,1,ni),[3,2,1]);
V12y = permute(V12(2,1,ni),[3,2,1]);
V12z = permute(V12(3,1,ni),[3,2,1]);
plot3(V12x,V12y,V12z,'^b');
% Ankle axis position in shank segment
% V13: virtual marker 1 of segment 3
V13 = Mprod_array3(repmat(([Segment(3).nV(1,1)*eye(3),...
    (1 + Segment(3).nV(2,1))*eye(3), ...
    - Segment(3).nV(2,1)*eye(3), ...
    Segment(3).nV(3,1)*eye(3)]),[1 1 n]),Segment(3).Q);
V13x = permute(V13(1,1,ni),[3,2,1]);
V13y = permute(V13(2,1,ni),[3,2,1]);
V13z = permute(V13(3,1,ni),[3,2,1]);
plot3(V13x,V13y,V13z,'^r'); % Should be supperimposed with V12
% Orientation subtalar axis direction in foot segment
% n12: virtual normal 1 of segment 2
n12 = Mprod_array3(repmat(([Segment(2).nn(1,1)*eye(3),...
    (Segment(2).nn(2,1))*eye(3), ...
    - Segment(2).nn(2,1)*eye(3), ...
    Segment(2).nn(3,1)*eye(3)]),[1 1 n]),Segment(2).Q);
n12x = permute(n12(1,1,ni),[3,2,1]);
n12y = permute(n12(2,1,ni),[3,2,1]);
n12z = permute(n12(3,1,ni),[3,2,1]);
quiver3(V12x,V12y,V12z,n12x,n12y,n12z,0.15,'--k');
% Orientation of ankle axis direction in shank segment
% n13: virtual normal 1 of segment 3
n13 = Mprod_array3(repmat(([Segment(3).nn(1,1)*eye(3),...
    (Segment(3).nn(2,1))*eye(3), ...
    - Segment(3).nn(2,1)*eye(3), ...
    Segment(3).nn(3,1)*eye(3)]),[1 1 n]),Segment(3).Q);
n13x = permute(n13(1,1,ni),[3,2,1]);
n13y = permute(n13(2,1,ni),[3,2,1]);
n13z = permute(n13(3,1,ni),[3,2,1]);
quiver3(V13x,V13y,V13z,n13x,n13y,n13z,0.15,'--k');

% Knee
% Point of Tibia / Femur Medial Plane
% V43: virtual marker 4 of segment 3
V43 = Mprod_array3(repmat(([Segment(3).nV(1,4)*eye(3),...
    (1 + Segment(3).nV(2,4))*eye(3), ...
    - Segment(3).nV(2,4)*eye(3), ...
    Segment(3).nV(3,4)*eye(3)]),[1 1 n]),Segment(3).Q);  
V43x = permute(V43(1,1,ni),[3,2,1]);
V43y = permute(V43(2,1,ni),[3,2,1]);
V43z = permute(V43(3,1,ni),[3,2,1]);
plot3(V43x,V43y,V43z,'^r');
% Center of Tibia / Femur Medial Sphere
% V15: virtual marker 1 of segment 5
V15 = Mprod_array3(repmat(([Segment(5).nV(1,1)*eye(3),...
    (1 + Segment(5).nV(2,1))*eye(3), ...
    - Segment(5).nV(2,1)*eye(3), ...
    Segment(5).nV(3,1)*eye(3)]),[1 1 n]),Segment(5).Q);  
V15x = permute(V15(1,1,ni),[3,2,1]);
V15y = permute(V15(2,1,ni),[3,2,1]);
V15z = permute(V15(3,1,ni),[3,2,1]);
plot3(V15x,V15y,V15z,'^c');
plot3([V43x,V15x]',[V43y,V15y]',[V43z,V15z]','--k');
% Point of Tibia / Femur Lateral Plane
% V53: virtual marker 5 of segment 3
V53 = Mprod_array3(repmat(([Segment(3).nV(1,5)*eye(3),...
    (1 + Segment(3).nV(2,5))*eye(3), ...
    - Segment(3).nV(2,5)*eye(3), ...
    Segment(3).nV(3,5)*eye(3)]),[1 1 n]),Segment(3).Q);  
V53x = permute(V53(1,1,ni),[3,2,1]);
V53y = permute(V53(2,1,ni),[3,2,1]);
V53z = permute(V53(3,1,ni),[3,2,1]);
plot3(V53x,V53y,V53z,'^r');
% Center of Tibia / Femur Lateral Sphere
% V25: virtual marker 2 of segment 5
V25 = Mprod_array3(repmat(([Segment(5).nV(1,2)*eye(3),...
    (1 + Segment(5).nV(2,2))*eye(3), ...
    - Segment(5).nV(2,2)*eye(3), ...
    Segment(5).nV(3,2)*eye(3)]),[1 1 n]),Segment(5).Q);  
V25x = permute(V25(1,1,ni),[3,2,1]);
V25y = permute(V25(2,1,ni),[3,2,1]);
V25z = permute(V25(3,1,ni),[3,2,1]);
plot3(V25x,V25y,V25z,'^c');
plot3([V53x,V25x]',[V53y,V25y]',[V53z,V25z]','--k');
% Insertion of AC Ligament
% V63 virtual marker 6 of segment 3
V63 = Mprod_array3(repmat(([Segment(3).nV(1,6)*eye(3),...
    (1 + Segment(3).nV(2,6))*eye(3), ...
    - Segment(3).nV(2,6)*eye(3), ...
    Segment(3).nV(3,6)*eye(3)]),[1 1 n]),Segment(3).Q);  
V63x = permute(V63(1,1,ni),[3,2,1]);
V63y = permute(V63(2,1,ni),[3,2,1]);
V63z = permute(V63(3,1,ni),[3,2,1]);
plot3(V63x,V63y,V63z,'^r');
% Origin of ACL Ligament
% V35: virtual marker 3 of segment 5
V35 = Mprod_array3(repmat(([Segment(5).nV(1,3)*eye(3),...
    (1 + Segment(5).nV(2,3))*eye(3), ...
    - Segment(5).nV(2,3)*eye(3), ...
    Segment(5).nV(3,3)*eye(3)]),[1 1 n]),Segment(5).Q);  
V35x = permute(V35(1,1,ni),[3,2,1]);
V35y = permute(V35(2,1,ni),[3,2,1]);
V35z = permute(V35(3,1,ni),[3,2,1]);
plot3(V35x,V35y,V35z,'^c');
plot3([V63x,V35x]',[V63y,V35y]',[V63z,V35z]','--k');
% Insertion of PC Ligament
% V73 virtual marker 7 of segment 3
V73 = Mprod_array3(repmat(([Segment(3).nV(1,7)*eye(3),...
    (1 + Segment(3).nV(2,7))*eye(3), ...
    - Segment(3).nV(2,7)*eye(3), ...
    Segment(3).nV(3,7)*eye(3)]),[1 1 n]),Segment(3).Q);  
V73x = permute(V73(1,1,ni),[3,2,1]);
V73y = permute(V73(2,1,ni),[3,2,1]);
V73z = permute(V73(3,1,ni),[3,2,1]);
plot3(V73x,V73y,V73z,'^r');
% Origin of PCL Ligament
% V45: virtual marker 4 of segment 5
V45 = Mprod_array3(repmat(([Segment(5).nV(1,4)*eye(3),...
    (1 + Segment(5).nV(2,4))*eye(3), ...
    - Segment(5).nV(2,4)*eye(3), ...
    Segment(5).nV(3,4)*eye(3)]),[1 1 n]),Segment(5).Q);  
V45x = permute(V45(1,1,ni),[3,2,1]);
V45y = permute(V45(2,1,ni),[3,2,1]);
V45z = permute(V45(3,1,ni),[3,2,1]);
plot3(V45x,V45y,V45z,'^c');
plot3([V73x,V45x]',[V73y,V45y]',[V73z,V45z]','--k');
% Insertion of MC Ligament
% V83 virtual marker 8 of segment 3
V83 = Mprod_array3(repmat(([Segment(3).nV(1,8)*eye(3),...
    (1 + Segment(3).nV(2,8))*eye(3), ...
    - Segment(3).nV(2,8)*eye(3), ...
    Segment(3).nV(3,8)*eye(3)]),[1 1 n]),Segment(3).Q);  
V83x = permute(V83(1,1,ni),[3,2,1]);
V83y = permute(V83(2,1,ni),[3,2,1]);
V83z = permute(V83(3,1,ni),[3,2,1]);
plot3(V83x,V83y,V83z,'^r');
% Origin of MCL Ligament
% V55: virtual marker 5 of segment 5
V55 = Mprod_array3(repmat(([Segment(5).nV(1,5)*eye(3),...
    (1 + Segment(5).nV(2,5))*eye(3), ...
    - Segment(5).nV(2,5)*eye(3), ...
    Segment(5).nV(3,5)*eye(3)]),[1 1 n]),Segment(5).Q);  
V55x = permute(V55(1,1,ni),[3,2,1]);
V55y = permute(V55(2,1,ni),[3,2,1]);
V55z = permute(V55(3,1,ni),[3,2,1]);
plot3(V55x,V55y,V55z,'^c');
plot3([V83x,V55x]',[V83y,V55y]',[V83z,V55z]','--k');
% Orientation of the medial tibial plateau
n23 = Mprod_array3(repmat(([Segment(3).nn(1,2)*eye(3),...
    (Segment(3).nn(2,2))*eye(3), ...
    - Segment(3).nn(2,2)*eye(3), ...
    Segment(3).nn(3,2)*eye(3)]),[1 1 n]),Segment(3).Q);
n23x = permute(n23(1,1,ni),[3,2,1]);
n23y = permute(n23(2,1,ni),[3,2,1]);
n23z = permute(n23(3,1,ni),[3,2,1]);
quiver3(V43x,V43y,V43z,n23x,n23y,n23z,0.15,'--k');
% Orientation of the lateral tibial plateau
n33 = Mprod_array3(repmat(([Segment(3).nn(1,3)*eye(3),...
    (Segment(3).nn(2,3))*eye(3), ...
    - Segment(3).nn(2,3)*eye(3), ...
    Segment(3).nn(3,3)*eye(3)]),[1 1 n]),Segment(3).Q);
n33x = permute(n33(1,1,ni),[3,2,1]);
n33y = permute(n33(2,1,ni),[3,2,1]);
n33z = permute(n33(3,1,ni),[3,2,1]);
quiver3(V53x,V53y,V53z,n33x,n33y,n33z,0.15,'--k');

% Patella
% Insertion of PT Ligament
% V93: virtual marker 9 of segment 3
V93 = Mprod_array3(repmat(([Segment(3).nV(1,9)*eye(3),...
    (1 + Segment(3).nV(2,9))*eye(3), ...
    - Segment(3).nV(2,9)*eye(3), ...
    Segment(3).nV(3,9)*eye(3)]),[1 1 n]),Segment(3).Q);  
V93x = permute(V93(1,1,ni),[3,2,1]);
V93y = permute(V93(2,1,ni),[3,2,1]);
V93z = permute(V93(3,1,ni),[3,2,1]);
plot3(V93x,V93y,V93z,'^r');
% Origin of PT Ligament = D4
plot3([V93x,D4x]',[V93y,D4y]',[V93z,D4z]','--k');
% Point of patella hinge axis (patella)
% V14: virtual marker 1 of segment 4
V14 = Mprod_array3(repmat(([Segment(4).nV(1,1)*eye(3),...
    (1 + Segment(4).nV(2,1))*eye(3), ...
    - Segment(4).nV(2,1)*eye(3), ...
    Segment(4).nV(3,1)*eye(3)]),[1 1 n]),Segment(4).Q);
V14x = permute(V14(1,1,ni),[3,2,1]);
V14y = permute(V14(2,1,ni),[3,2,1]);
V14z = permute(V14(3,1,ni),[3,2,1]);
plot3(V14x,V14y,V14z,'^g');
% Point of patella hinge axis (thigh)
% V65: virtual marker 6 of segment 5
V65 = Mprod_array3(repmat(([Segment(5).nV(1,6)*eye(3),...
    (1 + Segment(5).nV(2,6))*eye(3), ...
    - Segment(5).nV(2,6)*eye(3), ...
    Segment(5).nV(3,6)*eye(3)]),[1 1 n]),Segment(5).Q);  
V65x = permute(V65(1,1,ni),[3,2,1]);
V65y = permute(V65(2,1,ni),[3,2,1]);
V65z = permute(V65(3,1,ni),[3,2,1]);
plot3(V65x,V65y,V65z,'^c'); % Should be supperimposed with V14
% Orientation of patella hinge axis (patella)
% n14: virtual normal 1 of segment 4
n14 = Mprod_array3(repmat(([Segment(4).nn(1,1)*eye(3),...
    (Segment(4).nn(2,1))*eye(3), ...
    - Segment(4).nn(2,1)*eye(3), ...
    Segment(4).nn(3,1)*eye(3)]),[1 1 n]),Segment(4).Q);
n14x = permute(n14(1,1,ni),[3,2,1]);
n14y = permute(n14(2,1,ni),[3,2,1]);
n14z = permute(n14(3,1,ni),[3,2,1]);
quiver3(V14x,V14y,V14z,n14x,n14y,n14z,0.15,'--k');
% Orientation of patella hinge axis (thigh)
% n15: virtual normal 1 of segment 5
n15 = Mprod_array3(repmat(([Segment(5).nn(1,1)*eye(3),...
    (Segment(5).nn(2,1))*eye(3), ...
    - Segment(5).nn(2,1)*eye(3), ...
    Segment(5).nn(3,1)*eye(3)]),[1 1 n]),Segment(5).Q);
n15x = permute(n15(1,1,ni),[3,2,1]);
n15y = permute(n15(2,1,ni),[3,2,1]);
n15z = permute(n15(3,1,ni),[3,2,1]);
quiver3(V65x,V65y,V65z,n15x,n15y,n15z,0.15,'--k');


% Hip
% Hip Joint Center
% V16: virtual marker 1 of segment 6
V16 = Mprod_array3(repmat(([Segment(6).nV(1,1)*eye(3),...
    (1 + Segment(6).nV(2,1))*eye(3), ...
    - Segment(6).nV(2,1)*eye(3), ...
    Segment(6).nV(3,1)*eye(3)]),[1 1 n]),Segment(6).Q);  
V16x = permute(V16(1,1,ni),[3,2,1]);
V16y = permute(V16(2,1,ni),[3,2,1]);
V16z = permute(V16(3,1,ni),[3,2,1]);
plot3(V16x,V16y,V16z,'^m');


%% Muscles

% Gluteus maximus I
V75 = Mprod_array3(repmat([Segment(5).nV(1,7)*eye(3), ...
    (1 + Segment(5).nV(2,7))*eye(3), ...
    - Segment(5).nV(2,7)*eye(3), ...
    Segment(5).nV(3,7)*eye(3)],[1 1 n]),Segment(5).Q);
V75x = permute(V75(1,1,ni),[3,2,1]);
V75y = permute(V75(2,1,ni),[3,2,1]);
V75z = permute(V75(3,1,ni),[3,2,1]);
plot3(V75x,V75y,V75z,'sc');
V26 = Mprod_array3(repmat([Segment(6).nV(1,2)*eye(3), ...
    (1 + Segment(6).nV(2,2))*eye(3), ...
    - Segment(6).nV(2,2)*eye(3), ...
    Segment(6).nV(3,2)*eye(3)],[1 1 n]),Segment(6).Q);
V26x = permute(V26(1,1,ni),[3,2,1]);
V26y = permute(V26(2,1,ni),[3,2,1]);
V26z = permute(V26(3,1,ni),[3,2,1]);
plot3(V26x,V26y,V26z,'sm');
plot3([V75x,V26x]',[V75y,V26y]',[V75z,V26z]','--k');

% Gluteus maximus II
V85 = Mprod_array3(repmat([Segment(5).nV(1,8)*eye(3), ...
    (1 + Segment(5).nV(2,8))*eye(3), ...
    - Segment(5).nV(2,8)*eye(3), ...
    Segment(5).nV(3,8)*eye(3)],[1 1 n]),Segment(5).Q);
V85x = permute(V85(1,1,ni),[3,2,1]);
V85y = permute(V85(2,1,ni),[3,2,1]);
V85z = permute(V85(3,1,ni),[3,2,1]);
plot3(V85x,V85y,V85z,'sc');
V36 = Mprod_array3(repmat([Segment(6).nV(1,3)*eye(3), ...
    (1 + Segment(6).nV(2,3))*eye(3), ...
    - Segment(6).nV(2,3)*eye(3), ...
    Segment(6).nV(3,3)*eye(3)],[1 1 n]),Segment(6).Q);
V36x = permute(V36(1,1,ni),[3,2,1]);
V36y = permute(V36(2,1,ni),[3,2,1]);
V36z = permute(V36(3,1,ni),[3,2,1]);
plot3(V36x,V36y,V36z,'sm');
plot3([V85x,V36x]',[V85y,V36y]',[V85z,V36z]','--k')

% Gluteus maximus III
V95 = Mprod_array3(repmat([Segment(5).nV(1,9)*eye(3), ...
    (1 + Segment(5).nV(2,9))*eye(3), ...
    - Segment(5).nV(2,9)*eye(3), ...
    Segment(5).nV(3,9)*eye(3)],[1 1 n]),Segment(5).Q);
V95x = permute(V95(1,1,ni),[3,2,1]);
V95y = permute(V95(2,1,ni),[3,2,1]);
V95z = permute(V95(3,1,ni),[3,2,1]);
plot3(V95x,V95y,V95z,'sc');
V46 = Mprod_array3(repmat([Segment(6).nV(1,4)*eye(3), ...
    (1 + Segment(6).nV(2,4))*eye(3), ...
    - Segment(6).nV(2,4)*eye(3), ...
    Segment(6).nV(3,4)*eye(3)],[1 1 n]),Segment(6).Q);
V46x = permute(V46(1,1,ni),[3,2,1]);
V46y = permute(V46(2,1,ni),[3,2,1]);
V46z = permute(V46(3,1,ni),[3,2,1]);
plot3(V46x,V46y,V46z,'sm');
plot3([V95x,V46x]',[V95y,V46y]',[V95z,V46z]','--k')

% Gluteus medius I
V105 = Mprod_array3(repmat([Segment(5).nV(1,10)*eye(3), ...
    (1 + Segment(5).nV(2,10))*eye(3), ...
    - Segment(5).nV(2,10)*eye(3), ...
    Segment(5).nV(3,10)*eye(3)],[1 1 n]),Segment(5).Q);
V105x = permute(V105(1,1,ni),[3,2,1]);
V105y = permute(V105(2,1,ni),[3,2,1]);
V105z = permute(V105(3,1,ni),[3,2,1]);
plot3(V105x,V105y,V105z,'sc');
V56 = Mprod_array3(repmat([Segment(6).nV(1,5)*eye(3), ...
    (1 + Segment(6).nV(2,5))*eye(3), ...
    - Segment(6).nV(2,5)*eye(3), ...
    Segment(6).nV(3,5)*eye(3)],[1 1 n]),Segment(6).Q);
V56x = permute(V56(1,1,ni),[3,2,1]);
V56y = permute(V56(2,1,ni),[3,2,1]);
V56z = permute(V56(3,1,ni),[3,2,1]);
plot3(V56x,V56y,V56z,'sm');
plot3([V105x,V56x]',[V105y,V56y]',[V105z,V56z]','--k')

% Gluteus medius II
V115 = Mprod_array3(repmat([Segment(5).nV(1,11)*eye(3), ...
    (1 + Segment(5).nV(2,11))*eye(3), ...
    - Segment(5).nV(2,11)*eye(3), ...
    Segment(5).nV(3,11)*eye(3)],[1 1 n]),Segment(5).Q);
V115x = permute(V115(1,1,ni),[3,2,1]);
V115y = permute(V115(2,1,ni),[3,2,1]);
V115z = permute(V115(3,1,ni),[3,2,1]);
plot3(V115x,V115y,V115z,'sc');
V66 = Mprod_array3(repmat([Segment(6).nV(1,6)*eye(3), ...
    (1 + Segment(6).nV(2,6))*eye(3), ...
    - Segment(6).nV(2,6)*eye(3), ...
    Segment(6).nV(3,6)*eye(3)],[1 1 n]),Segment(6).Q);
V66x = permute(V66(1,1,ni),[3,2,1]);
V66y = permute(V66(2,1,ni),[3,2,1]);
V66z = permute(V66(3,1,ni),[3,2,1]);
plot3(V66x,V66y,V66z,'sm');
plot3([V115x,V66x]',[V115y,V66y]',[V115z,V66z]','--k')

% Gluteus medius III
V125 = Mprod_array3(repmat([Segment(5).nV(1,12)*eye(3), ...
    (1 + Segment(5).nV(2,12))*eye(3), ...
    - Segment(5).nV(2,12)*eye(3), ...
    Segment(5).nV(3,12)*eye(3)],[1 1 n]),Segment(5).Q);
V125x = permute(V125(1,1,ni),[3,2,1]);
V125y = permute(V125(2,1,ni),[3,2,1]);
V125z = permute(V125(3,1,ni),[3,2,1]);
plot3(V125x,V125y,V125z,'sc');
V76 = Mprod_array3(repmat([Segment(6).nV(1,7)*eye(3), ...
    (1 + Segment(6).nV(2,7))*eye(3), ...
    - Segment(6).nV(2,7)*eye(3), ...
    Segment(6).nV(3,7)*eye(3)],[1 1 n]),Segment(6).Q);
V76x = permute(V76(1,1,ni),[3,2,1]);
V76y = permute(V76(2,1,ni),[3,2,1]);
V76z = permute(V76(3,1,ni),[3,2,1]);
plot3(V76x,V76y,V76z,'sm');
plot3([V125x,V76x]',[V125y,V76y]',[V125z,V76z]','--k')

% Gluteus minimus I
V135 = Mprod_array3(repmat([Segment(5).nV(1,13)*eye(3), ...
    (1 + Segment(5).nV(2,13))*eye(3), ...
    - Segment(5).nV(2,13)*eye(3), ...
    Segment(5).nV(3,13)*eye(3)],[1 1 n]),Segment(5).Q);
V135x = permute(V135(1,1,ni),[3,2,1]);
V135y = permute(V135(2,1,ni),[3,2,1]);
V135z = permute(V135(3,1,ni),[3,2,1]);
plot3(V135x,V135y,V135z,'sc');
V86 = Mprod_array3(repmat([Segment(6).nV(1,8)*eye(3), ...
    (1 + Segment(6).nV(2,8))*eye(3), ...
    - Segment(6).nV(2,8)*eye(3), ...
    Segment(6).nV(3,8)*eye(3)],[1 1 n]),Segment(6).Q);
V86x = permute(V86(1,1,ni),[3,2,1]);
V86y = permute(V86(2,1,ni),[3,2,1]);
V86z = permute(V86(3,1,ni),[3,2,1]);
plot3(V86x,V86y,V86z,'sm');
plot3([V135x,V86x]',[V135y,V86y]',[V135z,V86z]','--k')

% Gluteus minimus II
V145 = Mprod_array3(repmat([Segment(5).nV(1,14)*eye(3), ...
    (1 + Segment(5).nV(2,14))*eye(3), ...
    - Segment(5).nV(2,14)*eye(3), ...
    Segment(5).nV(3,14)*eye(3)],[1 1 n]),Segment(5).Q);
V145x = permute(V145(1,1,ni),[3,2,1]);
V145y = permute(V145(2,1,ni),[3,2,1]);
V145z = permute(V145(3,1,ni),[3,2,1]);
plot3(V145x,V145y,V145z,'sc');
V96 = Mprod_array3(repmat([Segment(6).nV(1,9)*eye(3), ...
    (1 + Segment(6).nV(2,9))*eye(3), ...
    - Segment(6).nV(2,9)*eye(3), ...
    Segment(6).nV(3,9)*eye(3)],[1 1 n]),Segment(6).Q);
V96x = permute(V96(1,1,ni),[3,2,1]);
V96y = permute(V96(2,1,ni),[3,2,1]);
V96z = permute(V96(3,1,ni),[3,2,1]);
plot3(V96x,V96y,V96z,'sm');
plot3([V145x,V96x]',[V145y,V96y]',[V145z,V96z]','--k')

% Gluteus minimus III
V155 = Mprod_array3(repmat([Segment(5).nV(1,15)*eye(3), ...
    (1 + Segment(5).nV(2,15))*eye(3), ...
    - Segment(5).nV(2,15)*eye(3), ...
    Segment(5).nV(3,15)*eye(3)],[1 1 n]),Segment(5).Q);
V155x = permute(V155(1,1,ni),[3,2,1]);
V155y = permute(V155(2,1,ni),[3,2,1]);
V155z = permute(V155(3,1,ni),[3,2,1]);
plot3(V155x,V155y,V155z,'sc');
V106 = Mprod_array3(repmat([Segment(6).nV(1,10)*eye(3), ...
    (1 + Segment(6).nV(2,10))*eye(3), ...
    - Segment(6).nV(2,10)*eye(3), ...
    Segment(6).nV(3,10)*eye(3)],[1 1 n]),Segment(6).Q);
V106x = permute(V106(1,1,ni),[3,2,1]);
V106y = permute(V106(2,1,ni),[3,2,1]);
V106z = permute(V106(3,1,ni),[3,2,1]);
plot3(V106x,V106y,V106z,'sm');
plot3([V155x,V106x]',[V155y,V106y]',[V155z,V106z]','--k')

% Adductor longus
V165 = Mprod_array3(repmat([Segment(5).nV(1,16)*eye(3), ...
    (1 + Segment(5).nV(2,16))*eye(3), ...
    - Segment(5).nV(2,16)*eye(3), ...
    Segment(5).nV(3,16)*eye(3)],[1 1 n]),Segment(5).Q);
V165x = permute(V165(1,1,ni),[3,2,1]);
V165y = permute(V165(2,1,ni),[3,2,1]);
V165z = permute(V165(3,1,ni),[3,2,1]);
plot3(V165x,V165y,V165z,'sc');
V116 = Mprod_array3(repmat([Segment(6).nV(1,11)*eye(3), ...
    (1 + Segment(6).nV(2,11))*eye(3), ...
    - Segment(6).nV(2,11)*eye(3), ...
    Segment(6).nV(3,11)*eye(3)],[1 1 n]),Segment(6).Q);
V116x = permute(V116(1,1,ni),[3,2,1]);
V116y = permute(V116(2,1,ni),[3,2,1]);
V116z = permute(V116(3,1,ni),[3,2,1]);
plot3(V116x,V116y,V116z,'sm');
plot3([V165x,V116x]',[V165y,V116y]',[V165z,V116z]','--k')

% Adductor brevis
V175 = Mprod_array3(repmat([Segment(5).nV(1,17)*eye(3), ...
    (1 + Segment(5).nV(2,17))*eye(3), ...
    - Segment(5).nV(2,17)*eye(3), ...
    Segment(5).nV(3,17)*eye(3)],[1 1 n]),Segment(5).Q);
V175x = permute(V175(1,1,ni),[3,2,1]);
V175y = permute(V175(2,1,ni),[3,2,1]);
V175z = permute(V175(3,1,ni),[3,2,1]);
plot3(V175x,V175y,V175z,'sc');
V126 = Mprod_array3(repmat([Segment(6).nV(1,12)*eye(3), ...
    (1 + Segment(6).nV(2,12))*eye(3), ...
    - Segment(6).nV(2,12)*eye(3), ...
    Segment(6).nV(3,12)*eye(3)],[1 1 n]),Segment(6).Q);
V126x = permute(V126(1,1,ni),[3,2,1]);
V126y = permute(V126(2,1,ni),[3,2,1]);
V126z = permute(V126(3,1,ni),[3,2,1]);
plot3(V126x,V126y,V126z,'sm');
plot3([V175x,V126x]',[V175y,V126y]',[V175z,V126z]','--k')

% Adductor magnus I
V185 = Mprod_array3(repmat([Segment(5).nV(1,18)*eye(3), ...
    (1 + Segment(5).nV(2,18))*eye(3), ...
    - Segment(5).nV(2,18)*eye(3), ...
    Segment(5).nV(3,18)*eye(3)],[1 1 n]),Segment(5).Q);
V185x = permute(V185(1,1,ni),[3,2,1]);
V185y = permute(V185(2,1,ni),[3,2,1]);
V185z = permute(V185(3,1,ni),[3,2,1]);
plot3(V185x,V185y,V185z,'sc');
V136 = Mprod_array3(repmat([Segment(6).nV(1,13)*eye(3), ...
    (1 + Segment(6).nV(2,13))*eye(3), ...
    - Segment(6).nV(2,13)*eye(3), ...
    Segment(6).nV(3,13)*eye(3)],[1 1 n]),Segment(6).Q);
V136x = permute(V136(1,1,ni),[3,2,1]);
V136y = permute(V136(2,1,ni),[3,2,1]);
V136z = permute(V136(3,1,ni),[3,2,1]);
plot3(V136x,V136y,V136z,'sm');
plot3([V185x,V136x]',[V185y,V136y]',[V185z,V136z]','--k')

% Adductor magnus II
V195 = Mprod_array3(repmat([Segment(5).nV(1,19)*eye(3), ...
    (1 + Segment(5).nV(2,19))*eye(3), ...
    - Segment(5).nV(2,19)*eye(3), ...
    Segment(5).nV(3,19)*eye(3)],[1 1 n]),Segment(5).Q);
V195x = permute(V195(1,1,ni),[3,2,1]);
V195y = permute(V195(2,1,ni),[3,2,1]);
V195z = permute(V195(3,1,ni),[3,2,1]);
plot3(V195x,V195y,V195z,'sc');
V146 = Mprod_array3(repmat([Segment(6).nV(1,14)*eye(3), ...
    (1 + Segment(6).nV(2,14))*eye(3), ...
    - Segment(6).nV(2,14)*eye(3), ...
    Segment(6).nV(3,14)*eye(3)],[1 1 n]),Segment(6).Q);
V146x = permute(V146(1,1,ni),[3,2,1]);
V146y = permute(V146(2,1,ni),[3,2,1]);
V146z = permute(V146(3,1,ni),[3,2,1]);
plot3(V146x,V146y,V146z,'sm');
plot3([V195x,V146x]',[V195y,V146y]',[V195z,V146z]','--k')

% Adductor magnus III
V205 = Mprod_array3(repmat([Segment(5).nV(1,20)*eye(3), ...
    (1 + Segment(5).nV(2,20))*eye(3), ...
    - Segment(5).nV(2,20)*eye(3), ...
    Segment(5).nV(3,20)*eye(3)],[1 1 n]),Segment(5).Q);
V205x = permute(V205(1,1,ni),[3,2,1]);
V205y = permute(V205(2,1,ni),[3,2,1]);
V205z = permute(V205(3,1,ni),[3,2,1]);
plot3(V205x,V205y,V205z,'sc');
V156 = Mprod_array3(repmat([Segment(6).nV(1,15)*eye(3), ...
    (1 + Segment(6).nV(2,15))*eye(3), ...
    - Segment(6).nV(2,15)*eye(3), ...
    Segment(6).nV(3,15)*eye(3)],[1 1 n]),Segment(6).Q);
V156x = permute(V156(1,1,ni),[3,2,1]);
V156y = permute(V156(2,1,ni),[3,2,1]);
V156z = permute(V156(3,1,ni),[3,2,1]);
plot3(V156x,V156y,V156z,'sm');
plot3([V205x,V156x]',[V205y,V156y]',[V205z,V156z]','--k')

% Pectineus
V215 = Mprod_array3(repmat([Segment(5).nV(1,21)*eye(3), ...
    (1 + Segment(5).nV(2,21))*eye(3), ...
    - Segment(5).nV(2,21)*eye(3), ...
    Segment(5).nV(3,21)*eye(3)],[1 1 n]),Segment(5).Q);
V215x = permute(V215(1,1,ni),[3,2,1]);
V215y = permute(V215(2,1,ni),[3,2,1]);
V215z = permute(V215(3,1,ni),[3,2,1]);
plot3(V215x,V215y,V215z,'sc');
V166 = Mprod_array3(repmat([Segment(6).nV(1,16)*eye(3), ...
    (1 + Segment(6).nV(2,16))*eye(3), ...
    - Segment(6).nV(2,16)*eye(3), ...
    Segment(6).nV(3,16)*eye(3)],[1 1 n]),Segment(6).Q);
V166x = permute(V166(1,1,ni),[3,2,1]);
V166y = permute(V166(2,1,ni),[3,2,1]);
V166z = permute(V166(3,1,ni),[3,2,1]);
plot3(V166x,V166y,V166z,'sm');
plot3([V215x,V166x]',[V215y,V166y]',[V215z,V166z]','--k')

% Illiacus
V225 = Mprod_array3(repmat([Segment(5).nV(1,22)*eye(3), ...
    (1 + Segment(5).nV(2,22))*eye(3), ...
    - Segment(5).nV(2,22)*eye(3), ...
    Segment(5).nV(3,22)*eye(3)],[1 1 n]),Segment(5).Q);
V225x = permute(V225(1,1,ni),[3,2,1]);
V225y = permute(V225(2,1,ni),[3,2,1]);
V225z = permute(V225(3,1,ni),[3,2,1]);
plot3(V225x,V225y,V225z,'sc');
V176 = Mprod_array3(repmat([Segment(6).nV(1,17)*eye(3), ...
    (1 + Segment(6).nV(2,17))*eye(3), ...
    - Segment(6).nV(2,17)*eye(3), ...
    Segment(6).nV(3,17)*eye(3)],[1 1 n]),Segment(6).Q);
V176x = permute(V176(1,1,ni),[3,2,1]);
V176y = permute(V176(2,1,ni),[3,2,1]);
V176z = permute(V176(3,1,ni),[3,2,1]);
plot3(V176x,V176y,V176z,'sm');
plot3([V225x,V176x]',[V225y,V176y]',[V225z,V176z]','--k')

% Psoas
V235 = Mprod_array3(repmat([Segment(5).nV(1,23)*eye(3), ...
    (1 + Segment(5).nV(2,23))*eye(3), ...
    - Segment(5).nV(2,23)*eye(3), ...
    Segment(5).nV(3,23)*eye(3)],[1 1 n]),Segment(5).Q);
V235x = permute(V235(1,1,ni),[3,2,1]);
V235y = permute(V235(2,1,ni),[3,2,1]);
V235z = permute(V235(3,1,ni),[3,2,1]);
plot3(V235x,V235y,V235z,'sc');
V186 = Mprod_array3(repmat([Segment(6).nV(1,18)*eye(3), ...
    (1 + Segment(6).nV(2,18))*eye(3), ...
    - Segment(6).nV(2,18)*eye(3), ...
    Segment(6).nV(3,18)*eye(3)],[1 1 n]),Segment(6).Q);
V186x = permute(V186(1,1,ni),[3,2,1]);
V186y = permute(V186(2,1,ni),[3,2,1]);
V186z = permute(V186(3,1,ni),[3,2,1]);
plot3(V186x,V186y,V186z,'sm');
plot3([V235x,V186x]',[V235y,V186y]',[V235z,V186z]','--k')

% Quadratus femoris
V245 = Mprod_array3(repmat([Segment(5).nV(1,24)*eye(3), ...
    (1 + Segment(5).nV(2,24))*eye(3), ...
    - Segment(5).nV(2,24)*eye(3), ...
    Segment(5).nV(3,24)*eye(3)],[1 1 n]),Segment(5).Q);
V245x = permute(V245(1,1,ni),[3,2,1]);
V245y = permute(V245(2,1,ni),[3,2,1]);
V245z = permute(V245(3,1,ni),[3,2,1]);
plot3(V245x,V245y,V245z,'sc');
V196 = Mprod_array3(repmat([Segment(6).nV(1,19)*eye(3), ...
    (1 + Segment(6).nV(2,19))*eye(3), ...
    - Segment(6).nV(2,19)*eye(3), ...
    Segment(6).nV(3,19)*eye(3)],[1 1 n]),Segment(6).Q);
V196x = permute(V196(1,1,ni),[3,2,1]);
V196y = permute(V196(2,1,ni),[3,2,1]);
V196z = permute(V196(3,1,ni),[3,2,1]);
plot3(V196x,V196y,V196z,'sm');
plot3([V245x,V196x]',[V245y,V196y]',[V245z,V196z]','--k')

% Gemellus
V255 = Mprod_array3(repmat([Segment(5).nV(1,25)*eye(3), ...
    (1 + Segment(5).nV(2,25))*eye(3), ...
    - Segment(5).nV(2,25)*eye(3), ...
    Segment(5).nV(3,25)*eye(3)],[1 1 n]),Segment(5).Q);
V255x = permute(V255(1,1,ni),[3,2,1]);
V255y = permute(V255(2,1,ni),[3,2,1]);
V255z = permute(V255(3,1,ni),[3,2,1]);
plot3(V255x,V255y,V255z,'sc');
V206 = Mprod_array3(repmat([Segment(6).nV(1,20)*eye(3), ...
    (1 + Segment(6).nV(2,20))*eye(3), ...
    - Segment(6).nV(2,20)*eye(3), ...
    Segment(6).nV(3,20)*eye(3)],[1 1 n]),Segment(6).Q);
V206x = permute(V206(1,1,ni),[3,2,1]);
V206y = permute(V206(2,1,ni),[3,2,1]);
V206z = permute(V206(3,1,ni),[3,2,1]);
plot3(V206x,V206y,V206z,'sm');
plot3([V255x,V206x]',[V255y,V206y]',[V255z,V206z]','--k')

% Piriformis
V265 = Mprod_array3(repmat([Segment(5).nV(1,26)*eye(3), ...
    (1 + Segment(5).nV(2,26))*eye(3), ...
    - Segment(5).nV(2,26)*eye(3), ...
    Segment(5).nV(3,26)*eye(3)],[1 1 n]),Segment(5).Q);
V265x = permute(V265(1,1,ni),[3,2,1]);
V265y = permute(V265(2,1,ni),[3,2,1]);
V265z = permute(V265(3,1,ni),[3,2,1]);
plot3(V265x,V265y,V265z,'sc');
V216 = Mprod_array3(repmat([Segment(6).nV(1,21)*eye(3), ...
    (1 + Segment(6).nV(2,21))*eye(3), ...
    - Segment(6).nV(2,21)*eye(3), ...
    Segment(6).nV(3,21)*eye(3)],[1 1 n]),Segment(6).Q);
V216x = permute(V216(1,1,ni),[3,2,1]);
V216y = permute(V216(2,1,ni),[3,2,1]);
V216z = permute(V216(3,1,ni),[3,2,1]);
plot3(V216x,V216y,V216z,'sm');
plot3([V265x,V216x]',[V265y,V216y]',[V265z,V216z]','--k')

% Tensor fasciae latae
V275 = Mprod_array3(repmat([Segment(5).nV(1,27)*eye(3), ...
    (1 + Segment(5).nV(2,27))*eye(3), ...
    - Segment(5).nV(2,27)*eye(3), ...
    Segment(5).nV(3,27)*eye(3)],[1 1 n]),Segment(5).Q);
V275x = permute(V275(1,1,ni),[3,2,1]);
V275y = permute(V275(2,1,ni),[3,2,1]);
V275z = permute(V275(3,1,ni),[3,2,1]);
plot3(V275x,V275y,V275z,'sc');
V226 = Mprod_array3(repmat([Segment(6).nV(1,22)*eye(3), ...
    (1 + Segment(6).nV(2,22))*eye(3), ...
    - Segment(6).nV(2,22)*eye(3), ...
    Segment(6).nV(3,22)*eye(3)],[1 1 n]),Segment(6).Q);
V226x = permute(V226(1,1,ni),[3,2,1]);
V226y = permute(V226(2,1,ni),[3,2,1]);
V226z = permute(V226(3,1,ni),[3,2,1]);
plot3(V226x,V226y,V226z,'sm');
plot3([V275x,V226x]',[V275y,V226y]',[V275z,V226z]','--k')
V103 = Mprod_array3(repmat([Segment(3).nV(1,10)*eye(3), ...
    (1 + Segment(3).nV(2,10))*eye(3), ...
    - Segment(3).nV(2,10)*eye(3), ...
    Segment(3).nV(3,10)*eye(3)],[1 1 n]),Segment(3).Q);
V103x = permute(V103(1,1,ni),[3,2,1]);
V103y = permute(V103(2,1,ni),[3,2,1]);
V103z = permute(V103(3,1,ni),[3,2,1]);
plot3(V103x,V103y,V103z,'sr');
V285 = Mprod_array3(repmat([Segment(5).nV(1,28)*eye(3), ...
    (1 + Segment(5).nV(2,28))*eye(3), ...
    - Segment(5).nV(2,28)*eye(3), ...
    Segment(5).nV(3,28)*eye(3)],[1 1 n]),Segment(5).Q);
V285x = permute(V285(1,1,ni),[3,2,1]);
V285y = permute(V285(2,1,ni),[3,2,1]);
V285z = permute(V285(3,1,ni),[3,2,1]);
plot3(V285x,V285y,V285z,'sc');
plot3([V103x,V285x]',[V103y,V285y]',[V103z,V285z]','--k')

% Gracilis
V113 = Mprod_array3(repmat([Segment(3).nV(1,11)*eye(3), ...
    (1 + Segment(3).nV(2,11))*eye(3), ...
    - Segment(3).nV(2,11)*eye(3), ...
    Segment(3).nV(3,11)*eye(3)],[1 1 n]),Segment(3).Q);
V113x = permute(V113(1,1,ni),[3,2,1]);
V113y = permute(V113(2,1,ni),[3,2,1]);
V113z = permute(V113(3,1,ni),[3,2,1]);
plot3(V113x,V113y,V113z,'sr');
V236 = Mprod_array3(repmat([Segment(6).nV(1,23)*eye(3), ...
    (1 + Segment(6).nV(2,23))*eye(3), ...
    - Segment(6).nV(2,23)*eye(3), ...
    Segment(6).nV(3,23)*eye(3)],[1 1 n]),Segment(6).Q);
V236x = permute(V236(1,1,ni),[3,2,1]);
V236y = permute(V236(2,1,ni),[3,2,1]);
V236z = permute(V236(3,1,ni),[3,2,1]);
plot3(V236x,V236y,V236z,'sm');
plot3([V113x,V236x]',[V113y,V236y]',[V113z,V236z]','--k')

% Sartorius
V153 = Mprod_array3(repmat([Segment(3).nV(1,15)*eye(3), ...
    (1 + Segment(3).nV(2,15))*eye(3), ...
    - Segment(3).nV(2,15)*eye(3), ...
    Segment(3).nV(3,15)*eye(3)],[1 1 n]),Segment(3).Q);
V153x = permute(V153(1,1,ni),[3,2,1]);
V153y = permute(V153(2,1,ni),[3,2,1]);
V153z = permute(V153(3,1,ni),[3,2,1]);
plot3(V153x,V153y,V153z,'sr');
V295 = Mprod_array3(repmat([Segment(5).nV(1,29)*eye(3), ...
    (1 + Segment(5).nV(2,29))*eye(3), ...
    - Segment(5).nV(2,29)*eye(3), ...
    Segment(5).nV(3,29)*eye(3)],[1 1 n]),Segment(5).Q);
V295x = permute(V295(1,1,ni),[3,2,1]);
V295y = permute(V295(2,1,ni),[3,2,1]);
V295z = permute(V295(3,1,ni),[3,2,1]);
plot3(V295x,V295y,V295z,'sc');
plot3([V153x,V295x]',[V153y,V295y]',[V153z,V295z]','--k')
V276 = Mprod_array3(repmat([Segment(6).nV(1,27)*eye(3), ...
    (1 + Segment(6).nV(2,27))*eye(3), ...
    - Segment(6).nV(2,27)*eye(3), ...
    Segment(6).nV(3,27)*eye(3)],[1 1 n]),Segment(6).Q);
V276x = permute(V276(1,1,ni),[3,2,1]);
V276y = permute(V276(2,1,ni),[3,2,1]);
V276z = permute(V276(3,1,ni),[3,2,1]);
plot3(V276x,V276y,V276z,'sm');
plot3([V295x,V276x]',[V295y,V276y]',[V295z,V276z]','--k')

% Semimenbranosus
V123 = Mprod_array3(repmat([Segment(3).nV(1,12)*eye(3), ...
    (1 + Segment(3).nV(2,12))*eye(3), ...
    - Segment(3).nV(2,12)*eye(3), ...
    Segment(3).nV(3,12)*eye(3)],[1 1 n]),Segment(3).Q);
V123x = permute(V123(1,1,ni),[3,2,1]);
V123y = permute(V123(2,1,ni),[3,2,1]);
V123z = permute(V123(3,1,ni),[3,2,1]);
plot3(V123x,V123y,V123z,'sr');
V246 = Mprod_array3(repmat([Segment(6).nV(1,24)*eye(3), ...
    (1 + Segment(6).nV(2,24))*eye(3), ...
    - Segment(6).nV(2,24)*eye(3), ...
    Segment(6).nV(3,24)*eye(3)],[1 1 n]),Segment(6).Q);
V246x = permute(V246(1,1,ni),[3,2,1]);
V246y = permute(V246(2,1,ni),[3,2,1]);
V246z = permute(V246(3,1,ni),[3,2,1]);
plot3(V246x,V246y,V246z,'sm');
plot3([V123x,V246x]',[V123y,V246y]',[V123z,V246z]','--k')

% Semitendinus
V133 = Mprod_array3(repmat([Segment(3).nV(1,13)*eye(3), ...
    (1 + Segment(3).nV(2,13))*eye(3), ...
    - Segment(3).nV(2,13)*eye(3), ...
    Segment(3).nV(3,13)*eye(3)],[1 1 n]),Segment(3).Q);
V133x = permute(V133(1,1,ni),[3,2,1]);
V133y = permute(V133(2,1,ni),[3,2,1]);
V133z = permute(V133(3,1,ni),[3,2,1]);
plot3(V133x,V133y,V133z,'sr');
V256 = Mprod_array3(repmat([Segment(6).nV(1,25)*eye(3), ...
    (1 + Segment(6).nV(2,25))*eye(3), ...
    - Segment(6).nV(2,25)*eye(3), ...
    Segment(6).nV(3,25)*eye(3)],[1 1 n]),Segment(6).Q);
V256x = permute(V256(1,1,ni),[3,2,1]);
V256y = permute(V256(2,1,ni),[3,2,1]);
V256z = permute(V256(3,1,ni),[3,2,1]);
plot3(V256x,V256y,V256z,'sm');
plot3([V133x,V256x]',[V133y,V256y]',[V133z,V256z]','--k')

% Biceps femoris long head
V143 = Mprod_array3(repmat([Segment(3).nV(1,14)*eye(3), ...
    (1 + Segment(3).nV(2,14))*eye(3), ...
    - Segment(3).nV(2,14)*eye(3), ...
    Segment(3).nV(3,14)*eye(3)],[1 1 n]),Segment(3).Q);
V143x = permute(V143(1,1,ni),[3,2,1]);
V143y = permute(V143(2,1,ni),[3,2,1]);
V143z = permute(V143(3,1,ni),[3,2,1]);
plot3(V143x,V143y,V143z,'sr');
V266 = Mprod_array3(repmat([Segment(6).nV(1,26)*eye(3), ...
    (1 + Segment(6).nV(2,26))*eye(3), ...
    - Segment(6).nV(2,26)*eye(3), ...
    Segment(6).nV(3,26)*eye(3)],[1 1 n]),Segment(6).Q);
V266x = permute(V266(1,1,ni),[3,2,1]);
V266y = permute(V266(2,1,ni),[3,2,1]);
V266z = permute(V266(3,1,ni),[3,2,1]);
plot3(V266x,V266y,V266z,'sm');
plot3([V143x,V266x]',[V143y,V266y]',[V143z,V266z]','--k')

% Biceps femoris short head
V163 = Mprod_array3(repmat([Segment(3).nV(1,16)*eye(3), ...
    (1 + Segment(3).nV(2,16))*eye(3), ...
    - Segment(3).nV(2,16)*eye(3), ...
    Segment(3).nV(3,16)*eye(3)],[1 1 n]),Segment(3).Q);
V163x = permute(V163(1,1,ni),[3,2,1]);
V163y = permute(V163(2,1,ni),[3,2,1]);
V163z = permute(V163(3,1,ni),[3,2,1]);
plot3(V163x,V163y,V163z,'sr');
V305 = Mprod_array3(repmat([Segment(5).nV(1,30)*eye(3), ...
    (1 + Segment(5).nV(2,30))*eye(3), ...
    - Segment(5).nV(2,30)*eye(3), ...
    Segment(5).nV(3,30)*eye(3)],[1 1 n]),Segment(5).Q);
V305x = permute(V305(1,1,ni),[3,2,1]);
V305y = permute(V305(2,1,ni),[3,2,1]);
V305z = permute(V305(3,1,ni),[3,2,1]);
plot3(V305x,V305y,V305z,'sc');
plot3([V163x,V305x]',[V163y,V305y]',[V163z,V305z]','--k')

% Rectus femoris
V24 = Mprod_array3(repmat([Segment(4).nV(1,2)*eye(3), ...
    (1 + Segment(4).nV(2,2))*eye(3), ...
    - Segment(4).nV(2,2)*eye(3), ...
    Segment(4).nV(3,2)*eye(3)],[1 1 n]),Segment(4).Q);
V24x = permute(V24(1,1,ni),[3,2,1]);
V24y = permute(V24(2,1,ni),[3,2,1]);
V24z = permute(V24(3,1,ni),[3,2,1]);
plot3(V24x,V24y,V24z,'sg');
V286 = Mprod_array3(repmat([Segment(6).nV(1,28)*eye(3), ...
    (1 + Segment(6).nV(2,28))*eye(3), ...
    - Segment(6).nV(2,28)*eye(3), ...
    Segment(6).nV(3,28)*eye(3)],[1 1 n]),Segment(6).Q);
V286x = permute(V286(1,1,ni),[3,2,1]);
V286y = permute(V286(2,1,ni),[3,2,1]);
V286z = permute(V286(3,1,ni),[3,2,1]);
plot3(V286x,V286y,V286z,'sm');
plot3([V24x,V286x]',[V24y,V286y]',[V24z,V286z]','--k')


% Vastus medialis
V34 = Mprod_array3(repmat([Segment(4).nV(1,3)*eye(3), ...
    (1 + Segment(4).nV(2,3))*eye(3), ...
    - Segment(4).nV(2,3)*eye(3), ...
    Segment(4).nV(3,3)*eye(3)],[1 1 n]),Segment(4).Q);
V34x = permute(V34(1,1,ni),[3,2,1]);
V34y = permute(V34(2,1,ni),[3,2,1]);
V34z = permute(V34(3,1,ni),[3,2,1]);
plot3(V34x,V34y,V34z,'sg');
V325 = Mprod_array3(repmat([Segment(5).nV(1,32)*eye(3), ...
    (1 + Segment(5).nV(2,32))*eye(3), ...
    - Segment(5).nV(2,32)*eye(3), ...
    Segment(5).nV(3,32)*eye(3)],[1 1 n]),Segment(5).Q);
V325x = permute(V325(1,1,ni),[3,2,1]);
V325y = permute(V325(2,1,ni),[3,2,1]);
V325z = permute(V325(3,1,ni),[3,2,1]);
plot3(V325x,V325y,V325z,'sc');
plot3([V34x,V325x]',[V34y,V325y]',[V34z,V325z]','--k')

% Vastus intermedialis
V44 = Mprod_array3(repmat([Segment(4).nV(1,4)*eye(3), ...
    (1 + Segment(4).nV(2,4))*eye(3), ...
    - Segment(4).nV(2,4)*eye(3), ...
    Segment(4).nV(3,4)*eye(3)],[1 1 n]),Segment(4).Q);
V44x = permute(V44(1,1,ni),[3,2,1]);
V44y = permute(V44(2,1,ni),[3,2,1]);
V44z = permute(V44(3,1,ni),[3,2,1]);
plot3(V44x,V44y,V44z,'sg');
V355 = Mprod_array3(repmat([Segment(5).nV(1,35)*eye(3), ...
    (1 + Segment(5).nV(2,35))*eye(3), ...
    - Segment(5).nV(2,35)*eye(3), ...
    Segment(5).nV(3,35)*eye(3)],[1 1 n]),Segment(5).Q);
V355x = permute(V355(1,1,ni),[3,2,1]);
V355y = permute(V355(2,1,ni),[3,2,1]);
V355z = permute(V355(3,1,ni),[3,2,1]);
plot3(V355x,V355y,V355z,'sc');
plot3([V44x,V355x]',[V44y,V355y]',[V44z,V355z]','--k')

% Vastus lateralis
V54 = Mprod_array3(repmat([Segment(4).nV(1,5)*eye(3), ...
    (1 + Segment(4).nV(2,5))*eye(3), ...
    - Segment(4).nV(2,5)*eye(3), ...
    Segment(4).nV(3,5)*eye(3)],[1 1 n]),Segment(4).Q);
V54x = permute(V54(1,1,ni),[3,2,1]);
V54y = permute(V54(2,1,ni),[3,2,1]);
V54z = permute(V54(3,1,ni),[3,2,1]);
plot3(V54x,V54y,V54z,'sg');
V375 = Mprod_array3(repmat([Segment(5).nV(1,37)*eye(3), ...
    (1 + Segment(5).nV(2,37))*eye(3), ...
    - Segment(5).nV(2,37)*eye(3), ...
    Segment(5).nV(3,37)*eye(3)],[1 1 n]),Segment(5).Q);
V375x = permute(V375(1,1,ni),[3,2,1]);
V375y = permute(V375(2,1,ni),[3,2,1]);
V375z = permute(V375(3,1,ni),[3,2,1]);
plot3(V375x,V375y,V375z,'sc');
plot3([V54x,V375x]',[V54y,V375y]',[V54z,V375z]','--k')

% Gastrocnemius medialis
V42 = Mprod_array3(repmat([Segment(2).nV(1,4)*eye(3), ...
    (1 + Segment(2).nV(2,4))*eye(3), ...
    - Segment(2).nV(2,4)*eye(3), ...
    Segment(2).nV(3,4)*eye(3)],[1 1 n]),Segment(2).Q);
V42x = permute(V42(1,1,ni),[3,2,1]);
V42y = permute(V42(2,1,ni),[3,2,1]);
V42z = permute(V42(3,1,ni),[3,2,1]);
plot3(V42x,V42y,V42z,'sb');
V173 = Mprod_array3(repmat([Segment(3).nV(1,17)*eye(3), ...
    (1 + Segment(3).nV(2,17))*eye(3), ...
    - Segment(3).nV(2,17)*eye(3), ...
    Segment(3).nV(3,17)*eye(3)],[1 1 n]),Segment(3).Q);
V173x = permute(V173(1,1,ni),[3,2,1]);
V173y = permute(V173(2,1,ni),[3,2,1]);
V173z = permute(V173(3,1,ni),[3,2,1]);
plot3(V173x,V173y,V173z,'sr');
plot3([V42x,V173x]',[V42y,V173y]',[V42z,V173z]','--k')
V405 = Mprod_array3(repmat([Segment(5).nV(1,40)*eye(3), ...
    (1 + Segment(5).nV(2,40))*eye(3), ...
    - Segment(5).nV(2,40)*eye(3), ...
    Segment(5).nV(3,40)*eye(3)],[1 1 n]),Segment(5).Q);
V405x = permute(V405(1,1,ni),[3,2,1]);
V405y = permute(V405(2,1,ni),[3,2,1]);
V405z = permute(V405(3,1,ni),[3,2,1]);
plot3(V405x,V405y,V405z,'sc');
plot3([V173x,V405x]',[V173y,V405y]',[V173z,V405z]','--k')

% Gastrocnemius lateralis   
V52 = Mprod_array3(repmat([Segment(2).nV(1,5)*eye(3), ...
    (1 + Segment(2).nV(2,5))*eye(3), ...
    - Segment(2).nV(2,5)*eye(3), ...
    Segment(2).nV(3,5)*eye(3)],[1 1 n]),Segment(2).Q);
V52x = permute(V52(1,1,ni),[3,2,1]);
V52y = permute(V52(2,1,ni),[3,2,1]);
V52z = permute(V52(3,1,ni),[3,2,1]);
plot3(V52x,V52y,V52z,'sb');
V183 = Mprod_array3(repmat([Segment(3).nV(1,18)*eye(3), ...
    (1 + Segment(3).nV(2,18))*eye(3), ...
    - Segment(3).nV(2,18)*eye(3), ...
    Segment(3).nV(3,18)*eye(3)],[1 1 n]),Segment(3).Q);
V183x = permute(V183(1,1,ni),[3,2,1]);
V183y = permute(V183(2,1,ni),[3,2,1]);
V183z = permute(V183(3,1,ni),[3,2,1]);
plot3(V183x,V183y,V183z,'sr');
plot3([V52x,V183x]',[V52y,V183y]',[V52z,V183z]','--k')
V425 = Mprod_array3(repmat([Segment(5).nV(1,42)*eye(3), ...
    (1 + Segment(5).nV(2,42))*eye(3), ...
    - Segment(5).nV(2,42)*eye(3), ...
    Segment(5).nV(3,42)*eye(3)],[1 1 n]),Segment(5).Q);
V425x = permute(V425(1,1,ni),[3,2,1]);
V425y = permute(V425(2,1,ni),[3,2,1]);
V425z = permute(V425(3,1,ni),[3,2,1]);
plot3(V425x,V425y,V425z,'sc');
plot3([V183x,V425x]',[V183y,V425y]',[V183z,V425z]','--k')

% Soleus
V62 = Mprod_array3(repmat([Segment(2).nV(1,6)*eye(3), ...
    (1 + Segment(2).nV(2,6))*eye(3), ...
    - Segment(2).nV(2,6)*eye(3), ...
    Segment(2).nV(3,6)*eye(3)],[1 1 n]),Segment(2).Q);
V62x = permute(V62(1,1,ni),[3,2,1]);
V62y = permute(V62(2,1,ni),[3,2,1]);
V62z = permute(V62(3,1,ni),[3,2,1]);
plot3(V62x,V62y,V62z,'sb');
V193 = Mprod_array3(repmat([Segment(3).nV(1,19)*eye(3), ...
    (1 + Segment(3).nV(2,19))*eye(3), ...
    - Segment(3).nV(2,19)*eye(3), ...
    Segment(3).nV(3,19)*eye(3)],[1 1 n]),Segment(3).Q);
V193x = permute(V193(1,1,ni),[3,2,1]);
V193y = permute(V193(2,1,ni),[3,2,1]);
V193z = permute(V193(3,1,ni),[3,2,1]);
plot3(V193x,V193y,V193z,'sr');
plot3([V62x,V193x]',[V62y,V193y]',[V62z,V193z]','--k')

% Tibialis posterior
V72 = Mprod_array3(repmat([Segment(2).nV(1,7)*eye(3), ...
    (1 + Segment(2).nV(2,7))*eye(3), ...
    - Segment(2).nV(2,7)*eye(3), ...
    Segment(2).nV(3,7)*eye(3)],[1 1 n]),Segment(2).Q);
V72x = permute(V72(1,1,ni),[3,2,1]);
V72y = permute(V72(2,1,ni),[3,2,1]);
V72z = permute(V72(3,1,ni),[3,2,1]);
plot3(V72x,V72y,V72z,'sb');
V203 = Mprod_array3(repmat([Segment(3).nV(1,20)*eye(3), ...
    (1 + Segment(3).nV(2,20))*eye(3), ...
    - Segment(3).nV(2,20)*eye(3), ...
    Segment(3).nV(3,20)*eye(3)],[1 1 n]),Segment(3).Q);
V203x = permute(V203(1,1,ni),[3,2,1]);
V203y = permute(V203(2,1,ni),[3,2,1]);
V203z = permute(V203(3,1,ni),[3,2,1]);
plot3(V203x,V203y,V203z,'sr');
plot3([V72x,V203x]',[V72y,V203y]',[V72z,V203z]','--k')

% Tibialis anterior
V82 = Mprod_array3(repmat([Segment(2).nV(1,8)*eye(3), ...
    (1 + Segment(2).nV(2,8))*eye(3), ...
    - Segment(2).nV(2,8)*eye(3), ...
    Segment(2).nV(3,8)*eye(3)],[1 1 n]),Segment(2).Q);
V82x = permute(V82(1,1,ni),[3,2,1]);
V82y = permute(V82(2,1,ni),[3,2,1]);
V82z = permute(V82(3,1,ni),[3,2,1]);
plot3(V82x,V82y,V82z,'sb');
V213 = Mprod_array3(repmat([Segment(3).nV(1,21)*eye(3), ...
    (1 + Segment(3).nV(2,21))*eye(3), ...
    - Segment(3).nV(2,21)*eye(3), ...
    Segment(3).nV(3,21)*eye(3)],[1 1 n]),Segment(3).Q);
V213x = permute(V213(1,1,ni),[3,2,1]);
V213y = permute(V213(2,1,ni),[3,2,1]);
V213z = permute(V213(3,1,ni),[3,2,1]);
plot3(V213x,V213y,V213z,'sr');
plot3([V82x,V213x]',[V82y,V213y]',[V82z,V213z]','--k')

% Peroneus brevis
V92 = Mprod_array3(repmat([Segment(2).nV(1,9)*eye(3), ...
    (1 + Segment(2).nV(2,9))*eye(3), ...
    - Segment(2).nV(2,9)*eye(3), ...
    Segment(2).nV(3,9)*eye(3)],[1 1 n]),Segment(2).Q);
V92x = permute(V92(1,1,ni),[3,2,1]);
V92y = permute(V92(2,1,ni),[3,2,1]);
V92z = permute(V92(3,1,ni),[3,2,1]);
plot3(V92x,V92y,V92z,'sb');
V223 = Mprod_array3(repmat([Segment(3).nV(1,22)*eye(3), ...
    (1 + Segment(3).nV(2,22))*eye(3), ...
    - Segment(3).nV(2,22)*eye(3), ...
    Segment(3).nV(3,22)*eye(3)],[1 1 n]),Segment(3).Q);
V223x = permute(V223(1,1,ni),[3,2,1]);
V223y = permute(V223(2,1,ni),[3,2,1]);
V223z = permute(V223(3,1,ni),[3,2,1]);
plot3(V223x,V223y,V223z,'sr');
plot3([V92x,V223x]',[V92y,V223y]',[V92z,V223z]','--k')

% Peroneus longus
V102 = Mprod_array3(repmat([Segment(2).nV(1,10)*eye(3), ...
    (1 + Segment(2).nV(2,10))*eye(3), ...
    - Segment(2).nV(2,10)*eye(3), ...
    Segment(2).nV(3,10)*eye(3)],[1 1 n]),Segment(2).Q);
V102x = permute(V102(1,1,ni),[3,2,1]);
V102y = permute(V102(2,1,ni),[3,2,1]);
V102z = permute(V102(3,1,ni),[3,2,1]);
plot3(V102x,V102y,V102z,'sb');
V233 = Mprod_array3(repmat([Segment(3).nV(1,23)*eye(3), ...
    (1 + Segment(3).nV(2,23))*eye(3), ...
    - Segment(3).nV(2,23)*eye(3), ...
    Segment(3).nV(3,23)*eye(3)],[1 1 n]),Segment(3).Q);
V233x = permute(V233(1,1,ni),[3,2,1]);
V233y = permute(V233(2,1,ni),[3,2,1]);
V233z = permute(V233(3,1,ni),[3,2,1]);
plot3(V233x,V233y,V233z,'sr');
plot3([V102x,V233x]',[V102y,V233y]',[V102z,V233z]','--k')

% Peroneus tertius
V112 = Mprod_array3(repmat([Segment(2).nV(1,11)*eye(3), ...
    (1 + Segment(2).nV(2,11))*eye(3), ...
    - Segment(2).nV(2,11)*eye(3), ...
    Segment(2).nV(3,11)*eye(3)],[1 1 n]),Segment(2).Q);
V112x = permute(V112(1,1,ni),[3,2,1]);
V112y = permute(V112(2,1,ni),[3,2,1]);
V112z = permute(V112(3,1,ni),[3,2,1]);
plot3(V112x,V112y,V112z,'sb');
V243 = Mprod_array3(repmat([Segment(3).nV(1,24)*eye(3),...
    (1 + Segment(3).nV(2,24))*eye(3), ...
    - Segment(3).nV(2,24)*eye(3), ...
    Segment(3).nV(3,24)*eye(3)],[1 1 n]),Segment(3).Q);
V243x = permute(V243(1,1,ni),[3,2,1]);
V243y = permute(V243(2,1,ni),[3,2,1]);
V243z = permute(V243(3,1,ni),[3,2,1]);
plot3(V243x,V243y,V243z,'sr');
plot3([V112x,V243x]',[V112y,V243y]',[V112z,V243z]','--k')

% Extensor digitorum longus
V122 = Mprod_array3(repmat([Segment(2).nV(1,12)*eye(3), ...
    (1 + Segment(2).nV(2,12))*eye(3), ...
    - Segment(2).nV(2,12)*eye(3), ...
    Segment(2).nV(3,12)*eye(3)],[1 1 n]),Segment(2).Q);
V122x = permute(V122(1,1,ni),[3,2,1]);
V122y = permute(V122(2,1,ni),[3,2,1]);
V122z = permute(V122(3,1,ni),[3,2,1]);
plot3(V122x,V122y,V122z,'sb');
V253 = Mprod_array3(repmat([Segment(3).nV(1,25)*eye(3), ...
    (1 + Segment(3).nV(2,25))*eye(3), ...
    - Segment(3).nV(2,25)*eye(3), ...
    Segment(3).nV(3,25)*eye(3)],[1 1 n]),Segment(3).Q);
V253x = permute(V253(1,1,ni),[3,2,1]);
V253y = permute(V253(2,1,ni),[3,2,1]);
V253z = permute(V253(3,1,ni),[3,2,1]);
plot3(V253x,V253y,V253z,'sr');
plot3([V122x,V253x]',[V122y,V253y]',[V122z,V253z]','--k')

% Extensor hallucis longus
V132 = Mprod_array3(repmat([Segment(2).nV(1,13)*eye(3), ...
    (1 + Segment(2).nV(2,13))*eye(3), ...
    - Segment(2).nV(2,13)*eye(3), ...
    Segment(2).nV(3,13)*eye(3)],[1 1 n]),Segment(2).Q);
V132x = permute(V132(1,1,ni),[3,2,1]);
V132y = permute(V132(2,1,ni),[3,2,1]);
V132z = permute(V132(3,1,ni),[3,2,1]);
plot3(V132x,V132y,V132z,'sb');
V263 = Mprod_array3(repmat([Segment(3).nV(1,26)*eye(3), ...
    (1 + Segment(3).nV(2,26))*eye(3), ...
    - Segment(3).nV(2,26)*eye(3), ...
    Segment(3).nV(3,26)*eye(3)],[1 1 n]),Segment(3).Q);
V263x = permute(V263(1,1,ni),[3,2,1]);
V263y = permute(V263(2,1,ni),[3,2,1]);
V263z = permute(V263(3,1,ni),[3,2,1]);
plot3(V263x,V263y,V263z,'sr');
plot3([V132x,V263x]',[V132y,V263y]',[V132z,V263z]','--k')

% Flexor digitorum longus
V142 = Mprod_array3(repmat([Segment(2).nV(1,14)*eye(3), ...
    (1 + Segment(2).nV(2,14))*eye(3), ...
    - Segment(2).nV(2,14)*eye(3), ...
    Segment(2).nV(3,14)*eye(3)],[1 1 n]),Segment(2).Q);
V142x = permute(V142(1,1,ni),[3,2,1]);
V142y = permute(V142(2,1,ni),[3,2,1]);
V142z = permute(V142(3,1,ni),[3,2,1]);
plot3(V142x,V142y,V142z,'sb');
V273 = Mprod_array3(repmat([Segment(3).nV(1,27)*eye(3), ...
    (1 + Segment(3).nV(2,27))*eye(3), ...
    - Segment(3).nV(2,27)*eye(3), ...
    Segment(3).nV(3,27)*eye(3)],[1 1 n]),Segment(3).Q);
V273x = permute(V273(1,1,ni),[3,2,1]);
V273y = permute(V273(2,1,ni),[3,2,1]);
V273z = permute(V273(3,1,ni),[3,2,1]);
plot3(V273x,V273y,V273z,'sr');
plot3([V142x,V273x]',[V142y,V273y]',[V142z,V273z]','--k')

% Flexor hallucis longus
V152 = Mprod_array3(repmat([Segment(2).nV(1,15)*eye(3), ...
    (1 + Segment(2).nV(2,15))*eye(3), ...
    - Segment(2).nV(2,15)*eye(3), ...
    Segment(2).nV(3,15)*eye(3)],[1 1 n]),Segment(2).Q);
V152x = permute(V152(1,1,ni),[3,2,1]);
V152y = permute(V152(2,1,ni),[3,2,1]);
V152z = permute(V152(3,1,ni),[3,2,1]);
plot3(V152x,V152y,V152z,'sb');
V283 = Mprod_array3(repmat([Segment(3).nV(1,28)*eye(3), ...
    (1 + Segment(3).nV(2,28))*eye(3), ...
    - Segment(3).nV(2,28)*eye(3), ...
    Segment(3).nV(3,28)*eye(3)],[1 1 n]),Segment(3).Q);
V283x = permute(V283(1,1,ni),[3,2,1]);
V283y = permute(V283(2,1,ni),[3,2,1]);
V283z = permute(V283(3,1,ni),[3,2,1]);
plot3(V283x,V283y,V283z,'sr');
plot3([V152x,V283x]',[V152y,V283y]',[V152z,V283z]','--k')

