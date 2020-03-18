% MAIN PROGRAM
% Scale_Model.m
%__________________________________________________________________________
%
% PURPOSE
% Set generic model geometry, scale generic model
%
% SYNOPSIS
% (Segment,Model) = Scale_Model(Segment)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
% Subject.height (subject height in cm)
% Scaling (1: scaled by body height, 2: scaled by segments geometry)
%
% OUTPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% DESCRIPTION
% Set generic model geometry
% Scale generic model with respect to body height or segment geometry
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
%
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Florent Moissenet
% March 2020
%
%__________________________________________________________________________

function [Segment,Model] = Scale_Model(Segment,Model,subjectHeight,scaling)

if size(Segment,2) == 5 % Patella segment not defined
    
    Segment(4).L = mean(sqrt(sum((Segment(4).Q(4:6,1,:) - ...
        Segment(4).Q(7:9,1,:)).^2)),3); % Mean segment length
    if scaling == 1 % Scaled by body height
        Segment(4).rScale = subjectHeight/Model.geometry.height;
    elseif scaling == 2 % Scaled by segments geometry
        if (isfield(Segment,'W') && ~isempty(Segment(5).W)) % Scaled by segment width
            Segment(4).rScale = Segment(4).W/Model.geometry.Wm(5);
        else % Scaled by segment length
            Segment(4).rScale = Segment(4).L/Model.geometry.Lm(5);
        end
    end
    
else % Patella segment already defined
    
    % ---------------------------------------------------------------------
    % Subject geometry
    % ---------------------------------------------------------------------
    Segment(1).L = NaN; % No value for segment 1 (Forceplate)
    for i = 2:6 % From i = 2 (Foot) to i = 6 (Pelvis)
        if ~(isfield(Segment,'L') && ~isempty(Segment(i).L)) % Informed segment length
            Segment(i).L = mean(sqrt(sum((Segment(i).Q(4:6,1,:) - ...
                Segment(i).Q(7:9,1,:)).^2)),3); % Mean segment length
        end
    end

    % ---------------------------------------------------------------------
    % Homothety ratios
    % ---------------------------------------------------------------------
    for i = [3 5 2 4 6]
        if i == 4 || i == 6 % Use same scale ratio for patella (4) and pelvis (6) than for femur (5)
            Segment(i).rScale = Segment(5).rScale;
        elseif i == 2
            Segment(i).rScale = (Segment(3).rScale+Segment(5).rScale)/2; % Mean shank (3) and thigh (5) ratio applied on foot (2)
        else
            if scaling == 1 % Scaled by body height
                Segment(i).rScale = subjectHeight/Model.geometry.height;
            elseif scaling == 2 % Scaled by segments geometry
                if (isfield(Segment,'W') && ~isempty(Segment(i).W)) % Scaled by segment width
                    Segment(i).rScale = Segment(i).W/Model.geometry.Wm(i);
                else % Scaled by segment length
                    Segment(i).rScale = Segment(i).L/Model.geometry.Lm(i);
                end
            end
        end
    end
    
end