function gc = processing(grandChallenge,filename)

% -------------------------------------------------------------------------
% Load files
% -------------------------------------------------------------------------
cd('C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Multi-objective optimisation\data');

if grandChallenge == 1
    cd 'grand_challenge_1';
    footoff = 63;
elseif grandChallenge == 2
    cd 'grand_challenge_2';
    footoff = 64;
elseif grandChallenge == 3
    cd 'grand_challenge_3';
    footoff = 65;
elseif grandChallenge == 5
    cd 'grand_challenge_5';
    footoff = 68;
end

for i = 1:5
    gc.results(i) = load([filename{i},'_results.mat']);
end
gc.weight = gc.results(1).weight;
cd ..

% -------------------------------------------------------------------------
% Musculo-tendon forces
% -------------------------------------------------------------------------
% Phases
Phase1 = 1:10;      % loading responce
Phase2 = 11:30;     % midstance
Phase3 = 31:50;     % terminal stance
Phase4 = 51:footoff;     % pre-swing
Phase5 = footoff+1:73;     % initial swing
Phase6 = 74:87;     % midswing
Phase7 = 88:100;	% terminal Swing
phases = [0 10 30 50 footoff 73 87 100];

% Measurements
Emg = gc.results.Emg;
nEMG = 14;
EMG = zeros(nEMG,100);
MAX_EMG = zeros(nEMG,1);
MIN_EMG = zeros(nEMG,1);
n = length(Emg.gmax);
k = (1:n)';
k0 = (linspace(1,n,100))';

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

onoff_EMG = zeros(nEMG,7);
ratio_EMG = zeros(nEMG,1);

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
end

% Estimations (BW)
for i = 1:14
    F_est0 = [];
    F_est1 = [];
    F_est2 = [];
    for j = 1:5
        n = size(gc.results(j).Model0.Fm,3);
        timing = 1:100*n/100;
        if i == 1
            F_est0 = [F_est0 interpft(squeeze(sum(gc.results(j).Model0.Fm(1:3,:,timing))),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(sum(gc.results(j).Model1.Fm(1:3,:,timing))),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(sum(gc.results(j).Model2.Fm(1:3,:,timing))),101)/(9.81*gc.weight)];
        elseif i == 2
            F_est0 = [F_est0 interpft(squeeze(sum(gc.results(j).Model0.Fm(4:6,:,timing))),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(sum(gc.results(j).Model1.Fm(4:6,:,timing))),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(sum(gc.results(j).Model2.Fm(4:6,:,timing))),101)/(9.81*gc.weight)];
        elseif i == 3
            F_est0 = [F_est0 interpft(squeeze(sum(gc.results(j).Model0.Fm(12:14,:,timing))),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(sum(gc.results(j).Model1.Fm(12:14,:,timing))),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(sum(gc.results(j).Model2.Fm(12:14,:,timing))),101)/(9.81*gc.weight)];
        elseif i == 4
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(21,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(21,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(21,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 5
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(24,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(24,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(24,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 6
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(26,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(26,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(26,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 7
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(28,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(28,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(28,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 8
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(29,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(29,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(29,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 9
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(31,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(31,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(31,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 10
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(32,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(32,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(32,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 11
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(33,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(33,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(33,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 12
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(34,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(34,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(34,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 13
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(36,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(36,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(36,:,timing)),101)/(9.81*gc.weight)];
        elseif i == 14
            F_est0 = [F_est0 interpft(squeeze(gc.results(j).Model0.Fm(39,:,timing)),101)/(9.81*gc.weight)];
            F_est1 = [F_est1 interpft(squeeze(gc.results(j).Model1.Fm(39,:,timing)),101)/(9.81*gc.weight)];
            F_est2 = [F_est2 interpft(squeeze(gc.results(j).Model2.Fm(39,:,timing)),101)/(9.81*gc.weight)];
        end
    end
    gc.m(i).F_est0(:,:) = F_est0;
    gc.m(i).F_est1(:,:) = F_est1;
    gc.m(i).F_est2(:,:) = F_est2;   
end

% -------------------------------------------------------------------------
% Figures musculo-tendon forces
% -------------------------------------------------------------------------
for i = 4:11%1:14 only muscles crossing the knee
    figure('units','normalized','position',[0.3 0.3 .19 .25]);
    sub1 = subplot(2,1,1);
    hold on; box on; grid on;
    xlim([0,100]);
    if i == 1
       ylim([0,1]);
    elseif i == 2
       ylim([0,3]);
    elseif i == 3
       ylim([0,1]);
    elseif i == 4
       ylim([0,1.5]);
    elseif i == 5
       ylim([0,1.5]);
    elseif i == 6
       ylim([0,1.5]);
    elseif i == 7
       ylim([0,1.5]);
    elseif i == 8
       ylim([0,1.5]);
    elseif i == 9
       ylim([0,1.5]);
    elseif i == 10
       ylim([0,1.5]);
    elseif i == 11
       ylim([0,1.5]);
    elseif i == 12
       ylim([0,7]);
    elseif i == 13
       ylim([0,1]);
    elseif i == 14
       ylim([0,3]);
    end
    set(gca,'XTick',[0 20 40 60 80 100]);
    set(gca,'XTickLabel',{'' '' '' '' '' ''});
    set(gca,'YTick',[0 0.5 1 1.5]);
    line([footoff footoff],[0 10],'Color','black','Linestyle','--')
    corridor(mean(gc.m(i).F_est0(:,:),2),std(gc.m(i).F_est0(:,:),1,2),'red');
    corridor(mean(gc.m(i).F_est1(:,:),2),std(gc.m(i).F_est1(:,:),1,2),'green');
    corridor(mean(gc.m(i).F_est2(:,:),2),std(gc.m(i).F_est2(:,:),1,2),'blue');
    set(gca,'FontName','Times','FontSize',14);
    sub2 = subplot(2,1,2);
    hold on;
    xlim([0,100]);
    set(gca,'XTick',[0 20 40 60 80 100]);
    set(gca,'YTick',[]);
    for j = 1:7
        if onoff_EMG(i,j) == 1
            rectangle('Position',[phases(j),0.1,phases(j+1)-phases(j),0.01],'LineStyle','-','FaceColor',[0.7,0.7,0.7]);
        else
            rectangle('Position',[phases(j),0.1,phases(j+1)-phases(j),0.01],'LineStyle','-','FaceColor',[1,1,1]);
        end
    end
    set(gca,'FontName','Times','FontSize',14);
    set(sub1,'position',[0.19 0.20 0.8 0.7]);
    set(sub2,'position',[0.19 0.15 0.8 0.05]);
    saveas(gca, ['C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Multi-objective optimisation\results\grand_challenge_',num2str(grandChallenge),'_plot_',num2str(i)], 'fig');
    saveas(gca, ['C:\Users\florent.moissenet\Documents\Professionnel\publications\articles\1- en cours\Moissenet - Multi-objective optimisation\results\grand_challenge_',num2str(grandChallenge),'_plot_',num2str(i)], 'pdf');
end


% -------------------------------------------------------------------------
% Contact forces
% -------------------------------------------------------------------------
% Measurements (BW)
gc.Fmedial_mes = [];
gc.Flateral_mes = [];
gc.Ftotal_mes = [];
for i = 1:5
    n = length(gc.results(i).Force.KneeMedial);
    timing = 1:100*n/100;
    gc.Fmedial_mes = [gc.Fmedial_mes interpft(squeeze(gc.results(i).Force.KneeMedial(:,:,timing)),101)];
    gc.Flateral_mes = [gc.Flateral_mes interpft(squeeze(gc.results(i).Force.KneeLateral(:,:,timing)),101)];
    gc.Ftotal_mes = [gc.Ftotal_mes interpft(squeeze(gc.results(i).Force.KneeMedial(:,:,timing)),101)+interpft(squeeze(gc.results(i).Force.KneeLateral(:,:,timing)),101)];
end
gc.Fmedial_mes = (gc.Fmedial_mes)./(9.81*gc.weight);
gc.Flateral_mes = (gc.Flateral_mes)./(9.81*gc.weight);
gc.Ftotal_mes = (gc.Ftotal_mes)./(9.81*gc.weight);   

% Estimations (BW)
gc.Fmedial_est0 = [];
gc.Fmedial_est1 = [];
gc.Fmedial_est2 = [];
gc.Flateral_est0 = [];
gc.Flateral_est1 = [];
gc.Flateral_est2 = [];
gc.Ftotal_est0 = [];
gc.Ftotal_est1 = [];
gc.Ftotal_est2 = [];
for i = 1:5
    n = length(gc.results(i).Model1.Fc(4,:,:));
    timing = 1:100*n/100;
    gc.Fmedial_est0 = [gc.Fmedial_est0 interpft(squeeze(gc.results(i).Model0.Fc(4,:,timing)),101)];
    gc.Fmedial_est1 = [gc.Fmedial_est1 interpft(squeeze(gc.results(i).Model1.Fc(4,:,timing)),101)];
    gc.Fmedial_est2 = [gc.Fmedial_est2 interpft(squeeze(gc.results(i).Model2.Fc(4,:,timing)),101)];
    gc.Flateral_est0 = [gc.Flateral_est0 interpft(squeeze(gc.results(i).Model0.Fc(5,:,timing)),101)];
    gc.Flateral_est1 = [gc.Flateral_est1 interpft(squeeze(gc.results(i).Model1.Fc(5,:,timing)),101)];
    gc.Flateral_est2 = [gc.Flateral_est2 interpft(squeeze(gc.results(i).Model2.Fc(5,:,timing)),101)];
    gc.Ftotal_est0 = [gc.Ftotal_est0 interpft(squeeze(gc.results(i).Model0.Fc(4,:,timing)),101)+interpft(squeeze(gc.results(i).Model0.Fc(5,:,timing)),101)];
    gc.Ftotal_est1 = [gc.Ftotal_est1 interpft(squeeze(gc.results(i).Model1.Fc(4,:,timing)),101)+interpft(squeeze(gc.results(i).Model1.Fc(5,:,timing)),101)];
    gc.Ftotal_est2 = [gc.Ftotal_est2 interpft(squeeze(gc.results(i).Model2.Fc(4,:,timing)),101)+interpft(squeeze(gc.results(i).Model2.Fc(5,:,timing)),101)];
end
gc.Fmedial_est0 = (gc.Fmedial_est0)./(9.81*gc.weight);
gc.Fmedial_est1 = (gc.Fmedial_est1)./(9.81*gc.weight);
gc.Fmedial_est2 = (gc.Fmedial_est2)./(9.81*gc.weight);
gc.Flateral_est0 = (gc.Flateral_est0)./(9.81*gc.weight);
gc.Flateral_est1 = (gc.Flateral_est1)./(9.81*gc.weight);
gc.Flateral_est2 = (gc.Flateral_est2)./(9.81*gc.weight);
gc.Ftotal_est0 = (gc.Ftotal_est0)./(9.81*gc.weight);   
gc.Ftotal_est1 = (gc.Ftotal_est1)./(9.81*gc.weight);   
gc.Ftotal_est2 = (gc.Ftotal_est2)./(9.81*gc.weight);   

% Coefficient of concordance (%)
gc.concordance0 = [];
gc.concordance1 = [];
gc.concordance2 = [];
for i = 1:5
    gc.concordance0 = [gc.concordance0 gc.results(i).Model0.concordance_EMG];
    gc.concordance1 = [gc.concordance1 gc.results(i).Model1.concordance_EMG];
    gc.concordance2 = [gc.concordance2 gc.results(i).Model2.concordance_EMG];
end

% RMSE (BW) and R2 between estimated and measured tibiofemoral contact forces
for i = 1:5
    [gc.RMSE_Fmedial0(i) gc.R2_Fmedial0(i)] = goodnessFit(gc.Fmedial_mes(:,i),gc.Fmedial_est0(:,i));
    [gc.RMSE_Fmedial1(i) gc.R2_Fmedial1(i)] = goodnessFit(gc.Fmedial_mes(:,i),gc.Fmedial_est1(:,i));
    [gc.RMSE_Fmedial2(i) gc.R2_Fmedial2(i)] = goodnessFit(gc.Fmedial_mes(:,i),gc.Fmedial_est2(:,i));
    [gc.RMSE_Flateral0(i) gc.R2_Flateral0(i)] = goodnessFit(gc.Flateral_mes(:,i),gc.Flateral_est0(:,i));
    [gc.RMSE_Flateral1(i) gc.R2_Flateral1(i)] = goodnessFit(gc.Flateral_mes(:,i),gc.Flateral_est1(:,i));
    [gc.RMSE_Flateral2(i) gc.R2_Flateral2(i)] = goodnessFit(gc.Flateral_mes(:,i),gc.Flateral_est2(:,i));
    [gc.RMSE_Ftotal0(i) gc.R2_Ftotal0(i)] = goodnessFit(gc.Ftotal_mes(:,i),gc.Ftotal_est0(:,i));
    [gc.RMSE_Ftotal1(i) gc.R2_Ftotal1(i)] = goodnessFit(gc.Ftotal_mes(:,i),gc.Ftotal_est1(:,i));
    [gc.RMSE_Ftotal2(i) gc.R2_Ftotal2(i)] = goodnessFit(gc.Ftotal_mes(:,i),gc.Ftotal_est2(:,i));
end    

% -------------------------------------------------------------------------
% Figures contact forces
% -------------------------------------------------------------------------
figure;
subplot(3,1,1);
xlim([0 100]);
ylim([0 4.5]);
box on;
grid on;
hold on;
% title('Medial force');
corridor(mean(gc.Fmedial_mes,2),std(gc.Fmedial_mes,1,2),'black');
corridor(mean(gc.Fmedial_est0,2),std(gc.Fmedial_est0,1,2),'red');
corridor(mean(gc.Fmedial_est1,2),std(gc.Fmedial_est1,1,2),'green');
corridor(mean(gc.Fmedial_est2,2),std(gc.Fmedial_est2,1,2),'blue');
subplot(3,1,2);
xlim([0 100]);
ylim([0 4.5]);
box on;
grid on;
hold on;
% title('Lateral force');
corridor(mean(gc.Flateral_mes,2),std(gc.Flateral_mes,1,2),'black');
corridor(mean(gc.Flateral_est0,2),std(gc.Flateral_est0,1,2),'red');
corridor(mean(gc.Flateral_est1,2),std(gc.Flateral_est1,1,2),'green');
corridor(mean(gc.Flateral_est2,2),std(gc.Flateral_est2,1,2),'blue');
subplot(3,1,3);
xlim([0 100]);
ylim([0 7.5]);
box on;
grid on;
hold on;
% title('Total force');
corridor(mean(gc.Ftotal_mes,2),std(gc.Ftotal_mes,1,2),'black');
corridor(mean(gc.Ftotal_est0,2),std(gc.Ftotal_est0,1,2),'red');
corridor(mean(gc.Ftotal_est1,2),std(gc.Ftotal_est1,1,2),'green');
corridor(mean(gc.Ftotal_est2,2),std(gc.Ftotal_est2,1,2),'blue');


% -------------------------------------------------------------------------
% Ligament forces
% -------------------------------------------------------------------------
% Estimations (BW)
gc.Facl_est0 = [];
gc.Facl_est1 = [];
gc.Facl_est2 = [];
gc.Fpcl_est0 = [];
gc.Fpcl_est1 = [];
gc.Fpcl_est2 = [];
for i = 1:5
    n = length(gc.results(i).Model1.Fl(3,:,:));
    timing = 1:100*n/100;
    gc.Facl_est0 = [gc.Facl_est0 interpft(squeeze(gc.results(i).Model0.Fl(1,:,timing)),101)];
    gc.Facl_est1 = [gc.Facl_est1 interpft(squeeze(gc.results(i).Model1.Fl(1,:,timing)),101)];
    gc.Facl_est2 = [gc.Facl_est2 interpft(squeeze(gc.results(i).Model2.Fl(1,:,timing)),101)];
    gc.Fpcl_est0 = [gc.Fpcl_est0 interpft(squeeze(gc.results(i).Model0.Fl(2,:,timing)),101)];
    gc.Fpcl_est1 = [gc.Fpcl_est1 interpft(squeeze(gc.results(i).Model1.Fl(2,:,timing)),101)];
    gc.Fpcl_est2 = [gc.Fpcl_est2 interpft(squeeze(gc.results(i).Model2.Fl(2,:,timing)),101)];
end
gc.Facl_est0 = (gc.Facl_est0)./(9.81*gc.weight);
gc.Facl_est1 = (gc.Facl_est1)./(9.81*gc.weight);
gc.Facl_est2 = (gc.Facl_est2)./(9.81*gc.weight);
gc.Fpcl_est0 = (gc.Fpcl_est0)./(9.81*gc.weight);
gc.Fpcl_est1 = (gc.Fpcl_est1)./(9.81*gc.weight);
gc.Fpcl_est2 = (gc.Fpcl_est2)./(9.81*gc.weight);  

% -------------------------------------------------------------------------
% Figures ligament forces
% -------------------------------------------------------------------------
figure;
subplot(2,1,1);
xlim([0 100]);
% ylim([0 4.5]);
box on;
grid on;
hold on;
% title('ACL force');
corridor(mean(gc.Facl_est0,2),std(gc.Facl_est0,1,2),'red');
corridor(mean(gc.Facl_est1,2),std(gc.Facl_est1,1,2),'green');
corridor(mean(gc.Facl_est2,2),std(gc.Facl_est2,1,2),'blue');
subplot(2,1,2);
xlim([0 100]);
% ylim([0 4.5]);
box on;
grid on;
hold on;
% title('PCL force');
corridor(mean(gc.Fpcl_est0,2),std(gc.Fpcl_est0,1,2),'red');
corridor(mean(gc.Fpcl_est1,2),std(gc.Fpcl_est1,1,2),'green');
corridor(mean(gc.Fpcl_est2,2),std(gc.Fpcl_est2,1,2),'blue');