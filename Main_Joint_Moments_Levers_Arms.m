%% Coordinate partition

% Subtalar axis
Nn12 = [Segment(2).nn(1,1)*eye(3),...
    (Segment(2).nn(2,1))*eye(3), ...
    - Segment(2).nn(2,1)*eye(3), ...
    Segment(2).nn(3,1)*eye(3)];
n12 = Mprod_array3(repmat(Nn12,[1,1,n]),Segment(2).Q); 
% Ankle (talocrural) axis
Nn13 = [Segment(3).nn(1,1)*eye(3),...
    (Segment(3).nn(2,1))*eye(3), ...
    - Segment(3).nn(2,1)*eye(3), ...
    Segment(3).nn(3,1)*eye(3)];
n13 = Mprod_array3(repmat(Nn13,[1,1,n]),Segment(3).Q);

for k = 1:n
    
    % Coordinate partition
    Q = [Segment(2).Q(:,:,k); ... % Foot (12*1*n)
    Segment(3).Q(:,:,k); ... % Shank (12*1*n)
    Segment(4).Q(:,:,k); ... % Patella (12*1*n)
    Segment(5).Q(:,:,k); ... % Thigh (12*1*n)
    Segment(6).Q(:,:,k)]; % Pelvis (12*1*n)

    % Modify u and w axes to match the hinge axes
    Q(1:3,:) = n12(1:3,:,k); % u axis of foot
    Q(22:24,:) = n13(1:3,:,k); % w axis of shank
    
    % Numerical Jacobian of the map from Q to the DoF angles
    [partition(:,:,k)] = jacobianest(@Q2DoFs,Q); 
    % Parameter reduction 
    reduction(:,:,k) = inv([Model.K(:,:,k);partition(:,1:48,k)]); % Dynamic equations wihout Q parameters of pelvis
    Model.ZK(:,:,k) = reduction(:,end-5:end,k);
    
    % Constraint equations
    Aeq2(:,:,k) = Model.ZK(:,:,k)'*Model.Lever(:,:,k);
    Beq2(:,:,k) = Model.ZK(:,:,k)'*...
        (Model.G(:,:,k)*Model.d2Qdt2(:,:,k) - Model.P(:,:,k) - Model.R(:,:,k));
    
end


%% Joint moments
% b = Zk'*(GQ¨ - P - R)

figure
plot(squeeze(Beq2)')
ylabel('Moment (in N.m)')
xlabel('Sample instant of time')
legend({'DoF 1: Ankle F(+)/E(-)', ...
    'DoF 2: Ankle Ad(+)/Ab(-)', ...
    'DoF 3: Knee E(+)/F(-)', ...
    'DoF 4: Hip F(+)/E(-)', ...
    'DoF 5: Hip Ad(+)/Ab(-)', ...
    'DoF 6: Hip I(+)/E(-)R'})


%% Lever arms
% Zk'*L

% Ankle flexion-extension
T2 = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(3).Q)),Q2Tuv_array3(Segment(2).Q));
Euler2 = R2mobileZXY_array3(T2(1:3,1:3,:)); 
figure
plot(squeeze(Euler2(:,1,:))'*180/pi,squeeze(Aeq2(1,:,:))')

% % Knee flexion-extension
% T3 = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(5).Q)),Q2Tuv_array3(Segment(3).Q));
% Euler3 = R2mobileZXY_array3(T3(1:3,1:3,:)); 
% figure
% plot(-squeeze(Euler3(:,1,:))'*180/pi,-squeeze(Aeq2(3,:,:))')
% 
% % Hip flexion-extension
% T5 = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(6).Q)),Q2Tuv_array3(Segment(5).Q));
% Euler5 = R2mobileZXY_array3(T5(1:3,1:3,:)); 
% figure
% plot(squeeze(Euler5(:,1,:))'*180/pi,squeeze(Aeq2(4,:,:))')
% 
% % Hip adduction-abduction
% T5 = Mprod_array3(Tinv_array3(Q2Twu_array3(Segment(6).Q)),Q2Tuv_array3(Segment(5).Q));
% Euler5 = R2mobileZXY_array3(T5(1:3,1:3,:)); 
% figure
% plot(squeeze(Euler5(:,3,:))'*180/pi,squeeze(Aeq2(6,:,:))')

legend({'Gluteus maximus I', 'Gluteus maximus II', 'Gluteus maximus III', ...
    'Gluteus medius I', 'Gluteus medius II', 'Gluteus medius III', ...
    'Gluteus minimus I', 'Gluteus minimus II', 'Gluteus minimus III', ...
    'Adductor longus', 'Adductor brevis', ...
    'Adductor magnus I', 'Adductor magnus II', 'Adductor magnus III', ...
    'Pectineus', 'Illiacus', 'Psoas', ...
    'Quadratus femoris', 'Gemelli', ...
    'Piriformis','Tensor fasciae latae', ...
    'Gracilis', 'Sartorius', 'Semimembranosus', ...
    'Semitendinus', 'Biceps femoris long head', ...
    'Biceps femoris short head', 'Rectus femoris', ...
    'Vastus medialis', 'Vastus intermedialis', 'Vastus lateralis', ...
    'Gastrocnemius medialis', 'Gastrocnemius lateralis', ...
    'Soleus', 'Tibialis posterior', 'Tibialis anterior', ...
    'Peroneus brevis', 'Peroneus longus', 'Peroneus tertius', ...
    'Extensor digitorum longus', 'Extensor hallucis longus', ...
    'Flexor digitorum longus', 'Flexor hallucis longus'})
