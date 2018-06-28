function Model = Dickerson_new(Model,Emg,Geometry,n)

k = 1:n;
k0 = linspace(1,n,100);

% =========================================================================
% Define gait phases (Perry)
% =========================================================================
Phase1 = 1:10;      % loading responce
Phase2 = 11:30;     % midstance
Phase3 = 31:50;     % terminal stance
Phase4 = 51:60;     % pre-swing
Phase5 = 61:73;     % initial swing
Phase6 = 74:87;     % midswing
Phase7 = 88:100;	% terminal Swing

% =========================================================================
% Compute concordance to EMG measurements
% =========================================================================

% Extract EMG envelops (5 measures / side in this dataset)
% -------------------------------------------------------------------------
nEMG = 14;
EMG = zeros(nEMG,100);
MAX_EMG = zeros(nEMG,1);
MIN_EMG = zeros(nEMG,1);

EMG(1,:) = interp1(k,Emg.gmax(:)',k0,'pchip');        % gluteus_maximus
Emg.gmed(end-5:end) = Emg.gmed(end-11:end-6);
EMG(2,:) = interp1(k,Emg.gmed(:)',k0,'pchip');        % gluteus_medius
EMG(3,:) = interp1(k,Emg.addmagnus(:)',k0,'pchip');   % adductus_magnus
EMG(4,:) = interp1(k,Emg.tfl(:)',k0,'pchip');         % tensor_fascia_late
EMG(5,:) = interp1(k,Emg.semimem(:)',k0,'pchip');     % semimembranosus
EMG(6,:) = interp1(k,Emg.bifem(:)',k0,'pchip');       % biceps_femoris_long_head
EMG(7,:) = interp1(k,Emg.rf(:)',k0,'pchip');          % rectus_femoris
EMG(8,:) = interp1(k,Emg.vasmed(:)',k0,'pchip');      % vastus_medialis
EMG(9,:) = interp1(k,Emg.vaslat(:)',k0,'pchip');      % vastus_lateralis
EMG(10,:) = interp1(k,Emg.medgas(:)',k0,'pchip');     % gastrocnemius_medialis
EMG(11,:) = interp1(k,Emg.latgas(:)',k0,'pchip');     % gastrocnemius_lateralis
EMG(12,:) = interp1(k,Emg.soleus(:)',k0,'pchip');     % soleus
EMG(13,:) = interp1(k,Emg.tibant(:)',k0,'pchip');     % tibialis_anterior
EMG(14,:) = interp1(k,Emg.peronl(:)',k0,'pchip');     % peroneus_longus

MAX_EMG(1,:) = max(interp1(k,Emg.gmax(:)',k0,'pchip'));        % gluteus_maximus
MAX_EMG(2,:) = max(interp1(k,Emg.gmed(:)',k0,'pchip'));        % gluteus_medius
MAX_EMG(3,:) = max(interp1(k,Emg.addmagnus(:)',k0,'pchip'));   % adductus_magnus
MAX_EMG(4,:) = max(interp1(k,Emg.tfl(:)',k0,'pchip'));         % tensor_fascia_late
MAX_EMG(5,:) = max(interp1(k,Emg.semimem(:)',k0,'pchip'));     % semimembranosus
MAX_EMG(6,:) = max(interp1(k,Emg.bifem(:)',k0,'pchip'));       % biceps_femoris_long_head
MAX_EMG(7,:) = max(interp1(k,Emg.rf(:)',k0,'pchip'));          % rectus_femoris
MAX_EMG(8,:) = max(interp1(k,Emg.vasmed(:)',k0,'pchip'));      % vastus_medialis
MAX_EMG(9,:) = max(interp1(k,Emg.vaslat(:)',k0,'pchip'));      % vastus_lateralis
MAX_EMG(10,:) = max(interp1(k,Emg.medgas(:)',k0,'pchip'));     % gastrocnemius_medialis
MAX_EMG(11,:) = max(interp1(k,Emg.latgas(:)',k0,'pchip'));     % gastrocnemius_lateralis
MAX_EMG(12,:) = max(interp1(k,Emg.soleus(:)',k0,'pchip'));     % soleus
MAX_EMG(13,:) = max(interp1(k,Emg.tibant(:)',k0,'pchip'));     % tibialis_anterior
MAX_EMG(14,:) = max(interp1(k,Emg.peronl(:)',k0,'pchip'));     % peroneus_longus

MIN_EMG(1,:) = min(interp1(k,Emg.gmax(:)',k0,'pchip'));        % gluteus_maximus
MIN_EMG(2,:) = min(interp1(k,Emg.gmed(:)',k0,'pchip'));        % gluteus_medius
MIN_EMG(3,:) = min(interp1(k,Emg.addmagnus(:)',k0,'pchip'));   % adductus_magnus
MIN_EMG(4,:) = min(interp1(k,Emg.tfl(:)',k0,'pchip'));         % tensor_fascia_late
MIN_EMG(5,:) = min(interp1(k,Emg.semimem(:)',k0,'pchip'));     % semimembranosus
MIN_EMG(6,:) = min(interp1(k,Emg.bifem(:)',k0,'pchip'));       % biceps_femoris_long_head
MIN_EMG(7,:) = min(interp1(k,Emg.rf(:)',k0,'pchip'));          % rectus_femoris
MIN_EMG(8,:) = min(interp1(k,Emg.vasmed(:)',k0,'pchip'));      % vastus_medialis
MIN_EMG(9,:) = min(interp1(k,Emg.vaslat(:)',k0,'pchip'));      % vastus_lateralis
MIN_EMG(10,:) = min(interp1(k,Emg.medgas(:)',k0,'pchip'));     % gastrocnemius_medialis
MIN_EMG(11,:) = min(interp1(k,Emg.latgas(:)',k0,'pchip'));     % gastrocnemius_lateralis
MIN_EMG(12,:) = min(interp1(k,Emg.soleus(:)',k0,'pchip'));     % soleus
MIN_EMG(13,:) = min(interp1(k,Emg.tibant(:)',k0,'pchip'));     % tibialis_anterior
MIN_EMG(14,:) = min(interp1(k,Emg.peronl(:)',k0,'pchip'));     % peroneus_longus

% Extract estimated musculo-tendon forces associated to the previous EMG
% -------------------------------------------------------------------------
MTF = zeros(nEMG,100);
MAX_MTF = zeros(nEMG,1);
Model.X = squeeze(Model.X);

if strcmp(Geometry,'Delp')
    % Delp
    MTF(1,:) = interp1(k,sum(Model.X(1:3,:)),k0,'pchip');      % gluteus_maximus
    MTF(2,:) = interp1(k,sum(Model.X(4:6,:)),k0,'pchip');      % gluteus_medius
    MTF(3,:) = interp1(k,sum(Model.X(12:14,:)),k0,'pchip');    % adductus_magnus
    MTF(4,:) = interp1(k,Model.X(21,:),k0,'pchip');            % tensor_fascia_late
    MTF(5,:) = interp1(k,Model.X(24,:),k0,'pchip');            % semimembranosus
    MTF(6,:) = interp1(k,Model.X(26,:),k0,'pchip');            % biceps_femoris_long_head
    MTF(7,:) = interp1(k,Model.X(28,:),k0,'pchip');            % rectus_femoris
    MTF(8,:) = interp1(k,Model.X(29,:),k0,'pchip');            % vastus_medialis
    MTF(9,:) = interp1(k,Model.X(31,:),k0,'pchip');            % vastus_lateralis
    MTF(10,:) = interp1(k,Model.X(32,:),k0,'pchip');           % gastrocnemius_medialis
    MTF(11,:) = interp1(k,Model.X(33,:),k0,'pchip');           % gastrocnemius_lateralis
    MTF(12,:) = interp1(k,Model.X(34,:),k0,'pchip');           % soleus
    MTF(13,:) = interp1(k,Model.X(36,:),k0,'pchip');           % tibialis_anterior
    MTF(14,:) = interp1(k,Model.X(39,:),k0,'pchip');           % peroneus_longus

    MAX_MTF(1,:) = max(interp1(k,sum(Model.X(1:3,:)),k0,'pchip'));      % gluteus_maximus
    MAX_MTF(2,:) = max(interp1(k,sum(Model.X(4:6,:)),k0,'pchip'));      % gluteus_medius
    MAX_MTF(3,:) = max(interp1(k,sum(Model.X(12:14,:)),k0,'pchip'));    % adductus_magnus
    MAX_MTF(4,:) = max(interp1(k,Model.X(21,:),k0,'pchip'));            % tensor_fascia_late
    MAX_MTF(5,:) = max(interp1(k,Model.X(24,:),k0,'pchip'));            % semimembranosus
    MAX_MTF(6,:) = max(interp1(k,Model.X(26,:),k0,'pchip'));            % biceps_femoris_long_head
    MAX_MTF(7,:) = max(interp1(k,Model.X(28,:),k0,'pchip'));            % rectus_femoris
    MAX_MTF(8,:) = max(interp1(k,Model.X(29,:),k0,'pchip'));            % vastus_medialis
    MAX_MTF(9,:) = max(interp1(k,Model.X(31,:),k0,'pchip'));            % vastus_lateralis
    MAX_MTF(10,:) = max(interp1(k,Model.X(32,:),k0,'pchip'));           % gastrocnemius_medialis
    MAX_MTF(11,:) = max(interp1(k,Model.X(33,:),k0,'pchip'));           % gastrocnemius_lateralis
    MAX_MTF(12,:) = max(interp1(k,Model.X(34,:),k0,'pchip'));           % soleus
    MAX_MTF(13,:) = max(interp1(k,Model.X(36,:),k0,'pchip'));           % tibialis_anterior
    MAX_MTF(14,:) = max(interp1(k,Model.X(39,:),k0,'pchip'));           % peroneus_longus

elseif strcmp(Geometry,'KH')
    % Klein Horsman
    MTF(1,:) = interp1(k,sum(Model.X(38:49,:)),k0,'pchip');    % gluteus_maximus
    MTF(2,:) = interp1(k,sum(Model.X(50:61,:)),k0,'pchip');    % gluteus_medius
    MTF(3,:) = interp1(k,sum(Model.X(13:25,:)),k0,'pchip');    % adductus_magnus
    MTF(4,:) = interp1(k,sum(Model.X(102:103,:)),k0,'pchip');  % tensor_fascia_late
    MTF(5,:) = interp1(k,Model.X(94,:),k0,'pchip');            % semimembranosus
    MTF(6,:) = interp1(k,Model.X(26,:),k0,'pchip');            % biceps_femoris_long_head
    MTF(7,:) = interp1(k,sum(Model.X(90:91,:)),k0,'pchip');    % rectus_femoris
    MTF(8,:) = interp1(k,sum(Model.X(120:129,:)),k0,'pchip');  % vastus_medialis
    MTF(9,:) = interp1(k,sum(Model.X(112:119,:)),k0,'pchip');  % vastus_lateralis
    MTF(10,:) = interp1(k,Model.X(35,:),k0,'pchip');           % gastrocnemius_medialis
    MTF(11,:) = interp1(k,Model.X(34,:),k0,'pchip');           % gastrocnemius_lateralis
    MTF(12,:) = interp1(k,sum(Model.X(96:101,:)),k0,'pchip');  % soleus
    MTF(13,:) = interp1(k,Model.X(104,:),k0,'pchip');          % tibialis_anterior
    MTF(14,:) = interp1(k,Model.X(79,:),k0,'pchip');           % peroneus_longus

    MAX_MTF(1,:) = max(interp1(k,sum(Model.X(38:49,:)),k0,'pchip'));    % gluteus_maximus
    MAX_MTF(2,:) = max(interp1(k,sum(Model.X(50:61,:)),k0,'pchip'));    % gluteus_medius
    MAX_MTF(3,:) = max(interp1(k,sum(Model.X(13:25,:)),k0,'pchip'));    % adductus_magnus
    MAX_MTF(4,:) = max(interp1(k,sum(Model.X(102:103,:)),k0,'pchip'));  % tensor_fascia_late
    MAX_MTF(5,:) = max(interp1(k,Model.X(94,:),k0,'pchip'));            % semimembranosus
    MAX_MTF(6,:) = max(interp1(k,Model.X(26,:),k0,'pchip'));            % biceps_femoris_long_head
    MAX_MTF(7,:) = max(interp1(k,sum(Model.X(90:91,:)),k0,'pchip'));    % rectus_femoris
    MAX_MTF(8,:) = max(interp1(k,sum(Model.X(120:129,:)),k0,'pchip'));  % vastus_medialis
    MAX_MTF(9,:) = max(interp1(k,sum(Model.X(112:119,:)),k0,'pchip'));  % vastus_lateralis
    MAX_MTF(10,:) = max(interp1(k,Model.X(35,:),k0,'pchip'));           % gastrocnemius_medialis
    MAX_MTF(11,:) = max(interp1(k,Model.X(34,:),k0,'pchip'));           % gastrocnemius_lateralis
    MAX_MTF(12,:) = max(interp1(k,sum(Model.X(96:101,:)),k0,'pchip'));  % soleus
    MAX_MTF(13,:) = max(interp1(k,Model.X(104,:),k0,'pchip'));          % tibialis_anterior
    MAX_MTF(14,:) = max(interp1(k,Model.X(79,:),k0,'pchip'));           % peroneus_longus
end

% Build the on/off matrix 
% -------------------------------------------------------------------------
onoff_EMG = zeros(nEMG,7);
onoff_MTF = zeros(nEMG,7);
ratio_EMG = zeros(nEMG,1);
threshold_MTF = 0.10;

for m = 1:nEMG
    if m == 2
        threshold_EMG = 0.40; % gmed
    else
        threshold_EMG = 0.20;
    end
    % EMG matrix
    if mean(EMG(m,Phase1)) > threshold_EMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
        onoff_EMG(m,1) = 1;
    end
    if mean(EMG(m,Phase2)) > threshold_EMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
        onoff_EMG(m,2) = 1;
    end
    if mean(EMG(m,Phase3)) > threshold_EMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
        onoff_EMG(m,3) = 1;
    end
    if mean(EMG(m,Phase4)) > threshold_EMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
        onoff_EMG(m,4) = 1;
    end
    if mean(EMG(m,Phase5)) > threshold_EMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
        onoff_EMG(m,5) = 1;
    end
    if mean(EMG(m,Phase6)) > threshold_EMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
        onoff_EMG(m,6) = 1;
    end
    if mean(EMG(m,Phase7)) > threshold_EMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
        onoff_EMG(m,7) = 1;
    end
    % MTF matrix
    if mean(MTF(m,Phase1)) > threshold_MTF*MAX_MTF(m);
        onoff_MTF(m,1) = 1;
    end
    if mean(MTF(m,Phase2)) > threshold_MTF*MAX_MTF(m);
        onoff_MTF(m,2) = 1;
    end
    if mean(MTF(m,Phase3)) > threshold_MTF*MAX_MTF(m);
        onoff_MTF(m,3) = 1;
    end
    if mean(MTF(m,Phase4)) > threshold_MTF*MAX_MTF(m);
        onoff_MTF(m,4) = 1;
    end
    if mean(MTF(m,Phase5)) > threshold_MTF*MAX_MTF(m);
        onoff_MTF(m,5) = 1;
    end
    if mean(MTF(m,Phase6)) > threshold_MTF*MAX_MTF(m);
        onoff_MTF(m,6) = 1;
    end
    if mean(MTF(m,Phase7)) > threshold_MTF*MAX_MTF(m);
        onoff_MTF(m,7) = 1;
    end
end

% Compute the concordance Matrix
% 0 = concordance ; 1 = EMG:off MTF:on; -1 = EMF:on MTF:off
% -------------------------------------------------------------------------
concordance_matrix = onoff_MTF - onoff_EMG;

% Compute the concordance index (%)
% -------------------------------------------------------------------------
Model.concordance_EMG = length(find(concordance_matrix==0))*100/(7*nEMG);