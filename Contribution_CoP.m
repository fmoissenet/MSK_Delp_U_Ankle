Y = [zeros(1,1,61); ones(1,1,61); zeros(1,1,61)];
% frames 2 to 61

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

% Compute musculo-tendon contributions to CoP position
for i = 1:length(list)
    if i == 1
        eval(['F = Contribution.',list{i},'.F1R(:,:,1:61);']);
        eval(['M = Contribution.',list{i},'.M1R(:,:,1:61);']);
        F(find(abs(F)<1)) = 0;
        M(find(abs(M)<0.1)) = 0;
    else
        eval(['F = Contribution.muscleForce.',list{i},'.F1R(:,:,1:61);']);
        eval(['M = Contribution.muscleForce.',list{i},'.M1R(:,:,1:61);']);
        F(find(abs(F)<1)) = 0;
        M(find(abs(M)<0.1)) = 0;
    end
    % Avec seulement les deux premiers termes, on exprime ce qu'un muscle
    % pourrait induire comme modification de la trajectoire du CoP actuelle
    CoP = Mprod_array3(1/dot(F,F),cross(F,M)) - ...
        Mprod_array3(1/(dot(F,Y).*dot(F,F)), Mprod_array3(dot(F,M),cross(F,Y))) + ...
        Segment(1).Q(4:6,:,1:61);
    CoPx(:,i) = permute(CoP(1,:,:),[3,2,1]);
    CoPz(:,i) = permute(CoP(3,:,:),[3,2,1]);
end

% Measured CoP
CoPx(:,12) = permute(Segment(1).Q(4,:,1:61),[3,2,1]);
CoPz(:,12) = permute(Segment(1).Q(6,:,1:61),[3,2,1]);
CoPx(isinf(CoPx)) = NaN;
CoPz(isinf(CoPz)) = NaN;

for i = 1:11
    figure; hold on; axis equal;
    xlim([-0.10 0.10]);
    ylim([-0.03 0.23]);
    line([0 0],[-0.03 0.40],'Color','black','Linestyle','--');
    if i == 12
        title('CoP');
    else
        title(strrep(list{i},'O_',''));
    end
    
    % Compute the rotation matrix to realign the foot longitudinal axis
    vec = [(mean(Segment(2).rM(3,3,2:61),3)+mean(Segment(2).rM(3,2,2:61),3))/2-mean(Segment(2).rM(3,1,2:61),3); ...
        (mean(Segment(2).rM(1,2,2:61),3)+mean(Segment(2).rM(1,3,2:61),3))/2-mean(Segment(2).rM(1,1,2:61),3); ...
        0]/norm([(mean(Segment(2).rM(3,3,2:61),3)+mean(Segment(2).rM(3,2,2:61),3))/2-mean(Segment(2).rM(3,1,2:61),3); ...
        (mean(Segment(2).rM(1,2,2:61),3)+mean(Segment(2).rM(1,3,2:61),3))/2-mean(Segment(2).rM(1,1,2:61),3); ...
        0]);
    vec0 = [0;1;0];
    R = [dot(vec,vec0) -norm(cross(vec,vec0)) 0;...
        norm(cross(vec,vec0)) dot(vec,vec0) 0;...
        0 0 1];    
    
    % Plot contribution lines
    for j = 2:1:61
        Contribution.CoP(:,i,j) = R*[CoPz(j,i);CoPx(j,i);0]-R*[CoPz(2,12);CoPx(2,12);0];
        Contribution.CoP(:,12,j) = R*[CoPz(j,12);CoPx(j,12);0]-R*[CoPz(2,12);CoPx(2,12);0];
        if ~isnan(Contribution.CoP(1,i,j))
            if j<10 % loading response
                lcolor = 'red';
            elseif j>10 && j<30 % midstance
                lcolor = 'green';
            elseif j>30 && j<50 % terminal stance
                lcolor = 'blue';
            elseif j>50 && j<=61 % preswing
                lcolor = 'magenta';
            end
            line([Contribution.CoP(1,12,j) Contribution.CoP(1,i,j)],[Contribution.CoP(2,12,j) Contribution.CoP(2,i,j)],'Color',lcolor,'Linewidth',1);
        end
    end
    plot(squeeze(Contribution.CoP(1,12,2:61))',squeeze(Contribution.CoP(2,12,2:61))','Linestyle','-','Color','black','Linewidth',2);
    cd('C:\Users\florent.moissenet\Documents\Professionnel\publications\communications\2017\SB\Plots');
    saveas(gcf,[strrep(list{i},'O_',''),'.png']);   
end