function [Rn_est, Euler_est, Rn_err, Euler_err] = test_radio_navigation_UWB(Rn_ref, Euler_ref, dTheta_ref, Vb_ref, dt)
%  [Rn_, Vb_, Euler_, Rn_err, Vb_err, Euler_err] = test_navigation_UWB(Rn_ref, Vb_ref, Euler_ref, dThe_ref, dt)
%  Моделирует функционирование алгоритма системы навигации колесного робота
%  на базе локальной радионавигационной системы UWB с использованием
%  фильтра частиц (Particle Filter)
%
%   Входные аргументы:
%   Rn_ref - Опорные значения координат в навигационной СК
%   Vb_ref - Опорные значения скоростей в связанной СК
%   Euler_ref - Опорные значения углов ориентации
%   dThe_ref  - Опорные измерения датчиков угловой скорости 
%   dt - Шаг интегрированя
%
%
%   Выходные аргументы:
%   Rn_ - Оценки значений координат в навигационной СК
%   Euler_ - Оценки значений углов ориентации
%   Rn_err - Ошибки оценок значений координат в навигационной СК
%   Euler_err - Ошибки оценок значений углов ориентации

%% Инициализация анимации
%% 
animation_rate = 10;
figure;
ax = gca;
ax.FontSize = 14;
view(2);
axis equal;
hold on; grid on;
set(gca,'Xlim',[-60 60]);
set(gca,'Ylim',[-60 60]);
set(gca,'Zlim',[-5 5]);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
set(gcf,'renderer','opengl');
xlabel('X, м', 'FontSize', 16);
ylabel('Y, м', 'FontSize', 16);
tail_ref = plot(0,0,'k.-');
tail_est = plot(0,0,'b.-');
tail_length = 1800;

%% Оси опорного трехгранника
at1 = [0; 0; 2.5]*10;
at2 = [0; 2.5; 0]*10;
at3 = [2.5; 0; 0]*10;
st =  [0; 0; 0];
pat1 = plot3([st(1) at1(1)],[st(2) at1(2)],...
    [st(3) at1(3)],'b--','linewidth', 2);
pat2 = plot3([st(1) at2(1)],[st(2) at2(2)],...
    [st(3) at2(3)],'g--','linewidth', 2);
pat3 = plot3([st(1) at3(1)],[st(2) at3(2)],...
    [st(3) at3(3)],'r--','linewidth', 2);

%% Оси оценки
ai1 = [0; 0; 1.5;]*10;
ai2 = [0; 1.5; 0;]*10;
ai3 = [1.5; 0; 0;]*10;
sz =  [0; 0; 0;];
pai1 = plot3([sz(1) ai1(1)],[sz(2) ai1(2)],...
    [sz(3) ai1(3)],'b-','linewidth', 2);
pai2 = plot3([sz(1) ai2(1)],[sz(2) ai2(2)],...
    [sz(3) ai2(3)],'g-','linewidth', 2);
pai3 = plot3([sz(1) ai3(1)],[sz(2) ai3(2)],...
    [sz(3) ai3(3)],'r-','linewidth', 2);

%% 
Vert = [ 0 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 0 0 1 ; 1 0 1 ; 1 1 1 ; 0 1 1];
Vert(:,1) = (Vert(:,1)-0.5)*10;
Vert(:,2) = (Vert(:,2)-0.5)*10;
Vert(:,3) = (Vert(:,3)-0.5)*10;
Faces = [ 1 2 6 5 ; 2 3 7 6 ; 3 4 8 7 ; 4 1 5 8 ; 1 2 3 4 ; 5 6 7 8 ];
ptch.Vertices = Vert;
ptch.Faces = Faces;
ptch.FaceVertexCData = gray(6);
ptch.FaceColor = 'flat';
ptch.FaceVertexAlphaData = 0.2;
ptch.FaceAlpha = 'flat' ; 
patch_handle = patch(ptch);

handle_particles = plot(nan, nan, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

%% Базовые станции UWB
Anchors(1,:) = [ 50;  50;  0];
Anchors(2,:) = [ 50; -50;  0];
Anchors(3,:) = [-50;  50;  0];
Anchors(4,:) = [-50; -50;  0];

plot3(Anchors(1,1),Anchors(1,2),Anchors(1,3),'ko','MarkerSize', 15,'MarkerfaceColor','y',...
    'linewidth',2);
plot3(Anchors(2,1),Anchors(2,2),Anchors(2,3),'ko','MarkerSize', 15,'MarkerfaceColor','y',...
    'linewidth',2);
plot3(Anchors(3,1),Anchors(3,2),Anchors(3,3),'ko','MarkerSize', 15,'MarkerfaceColor','y',...
    'linewidth',2);
plot3(Anchors(4,1),Anchors(4,2),Anchors(4,3),'ko','MarkerSize', 15,'MarkerfaceColor','y',...
    'linewidth',2);

%% Анимация измерений UWB
pf1z = plot3([0 Anchors(1,1)],[0 Anchors(1,2)],...
    [0 Anchors(1,3)],'k-.','linewidth', 1);
pf2z = plot3([0 Anchors(2,1)],[0 Anchors(2,2)],...
    [0 Anchors(2,3)],'k-.','linewidth', 1);
pf3z = plot3([0 Anchors(3,1)],[0 Anchors(3,2)],...
    [0 Anchors(3,3)],'k-.','linewidth', 1);
pf4z = plot3([0 Anchors(4,1)],[0 Anchors(4,2)],...
    [0 Anchors(4,3)],'k-.','linewidth', 1);


%% Параметры UWB
uwb_noise = 5e-1;
uwb_rate  = 100;

%% Инициализация фильтра частиц
initPos = Rn_ref(1, 1:2)' + [5; 5];
initAng = Euler_ref(1, 1) + 0.1;
fltr = pf_init(initPos, initAng);

%% Логи 
Nsim = size(Rn_ref, 1);
Rn_est = zeros(Nsim, 3);
Euler_est = zeros(Nsim, 3);
Rn_err = zeros(Nsim, 3);
Euler_err = zeros(Nsim, 3);

%% Цикл моделирования
for t=1:Nsim

    % Управления
    dtheta = dTheta_ref(t, 1);
    upsilon = Vb_ref(t, 1);

    % Измерения UWB 
    pos_uwb = Rn_ref(t,:)';
    Ranges = [norm(pos_uwb-Anchors(1,:)');
        norm(pos_uwb-Anchors(2,:)');
        norm(pos_uwb-Anchors(3,:)');
        norm(pos_uwb-Anchors(4,:)')];
    Ranges = Ranges+randn(4, 1) * uwb_noise;
    if (mod(t, uwb_rate) == 0)
        uwb_update = true;
    else
        uwb_update = false;
    end

    %  Алгоритм навигации
    [fltr, robot_state] = radio_navigation_system_UWB(fltr, upsilon, dtheta, Ranges, Anchors, uwb_update, dt);

    % Логи
    Rn_est(t, :) = [robot_state(1:2, 1); 0];
    Euler_est(t, :) = [robot_state(3, 1); 0; 0];

    % Ошибки
    Cref = angle_to_dcm(Euler_ref(t, 1), 0, 0);
    C = angle_to_dcm(robot_state(3, 1), 0, 0);
    Cerr = Cref * C';
    [Euler_err(t, 1), Euler_err(t, 2), Euler_err(t, 3)] = dcm_to_angle(Cerr);
    Rn_err(t, :) = Rn_ref(t,:) - Rn_est(t, :);
    
    % Анимация
    if (mod(t, animation_rate) == 0)

        set(handle_particles, 'XData', fltr.p(1, :), 'YData', fltr.p(2, :));

        Cbnref = angle_to_dcm(Euler_ref(t, 1),Euler_ref(t, 2),Euler_ref(t, 3));
        a1_ = Cbnref'*at1+Rn_ref(t, :)';
        a2_ = Cbnref'*at2+Rn_ref(t, :)';
        a3_ = Cbnref'*at3+Rn_ref(t, :)';
        s_ = Rn_ref(t, :);
        set(pat1,'Xdata',[s_(1) a1_(1)]);
        set(pat1,'Ydata',[s_(2) a1_(2)]);
        set(pat1,'Zdata',[s_(3) a1_(3)]);
        set(pat2,'Xdata',[s_(1) a2_(1)]);
        set(pat2,'Ydata',[s_(2) a2_(2)]);
        set(pat2,'Zdata',[s_(3) a2_(3)]);
        set(pat3,'Xdata',[s_(1) a3_(1)]);
        set(pat3,'Ydata',[s_(2) a3_(2)]);
        set(pat3,'Zdata',[s_(3) a3_(3)]);

        Cbn = angle_to_dcm(robot_state(3, 1), 0, 0)';
        a1_ = Cbn * ai1 + Rn_est(t, :)';
        a2_ = Cbn * ai2 + Rn_est(t, :)';
        a3_ = Cbn * ai3 + Rn_est(t, :)';
        s_ = Rn_est(t, :)';
        set(pai1,'Xdata',[s_(1) a1_(1)]);
        set(pai1,'Ydata',[s_(2) a1_(2)]);
        set(pai1,'Zdata',[s_(3) a1_(3)]);
        set(pai2,'Xdata',[s_(1) a2_(1)]);
        set(pai2,'Ydata',[s_(2) a2_(2)]);
        set(pai2,'Zdata',[s_(3) a2_(3)]);
        set(pai3,'Xdata',[s_(1) a3_(1)]);
        set(pai3,'Ydata',[s_(2) a3_(2)]);
        set(pai3,'Zdata',[s_(3) a3_(3)]);

        Vert_ = Vert;
        for j=1:size(Vert,1)
            Vert_(j,:) = (Vert(j,:)*Cbn') + Rn_est(t, :);
        end
        set(patch_handle,'Vertices',Vert_);

        if (uwb_update == 1)
            set(pf1z,'Xdata',[Rn_ref(t, 1) Anchors(1,1)]);
            set(pf1z,'Ydata',[Rn_ref(t, 2) Anchors(1,2)]);
            set(pf1z,'Zdata',[Rn_ref(t, 3) Anchors(1,3)]);

            set(pf2z,'Xdata',[Rn_ref(t, 1) Anchors(2,1)]);
            set(pf2z,'Ydata',[Rn_ref(t, 2) Anchors(2,2)]);
            set(pf2z,'Zdata',[Rn_ref(t, 3) Anchors(2,3)]);

            set(pf3z,'Xdata',[Rn_ref(t, 1) Anchors(3,1)]);
            set(pf3z,'Ydata',[Rn_ref(t, 2) Anchors(3,2)]);
            set(pf3z,'Zdata',[Rn_ref(t, 3) Anchors(3,3)]);

            set(pf4z,'Xdata',[Rn_ref(t, 1) Anchors(4,1)]);
            set(pf4z,'Ydata',[Rn_ref(t, 2) Anchors(4,2)]);
            set(pf4z,'Zdata',[Rn_ref(t, 3) Anchors(4,3)]);
        else
            set(pf1z,'Xdata',[NaN NaN]);
            set(pf1z,'Ydata',[NaN NaN]);
            set(pf1z,'Zdata',[NaN NaN]);

            set(pf2z,'Xdata',[NaN NaN]);
            set(pf2z,'Ydata',[NaN NaN]);
            set(pf2z,'Zdata',[NaN NaN]);

            set(pf3z,'Xdata',[NaN NaN]);
            set(pf3z,'Ydata',[NaN NaN]);
            set(pf3z,'Zdata',[NaN NaN]);

            set(pf4z,'Xdata',[NaN NaN]);
            set(pf4z,'Ydata',[NaN NaN]);
            set(pf4z,'Zdata',[NaN NaN]);
        end    

        if (t > tail_length)
            set(tail_ref,...
                'Xdata',Rn_ref(t-tail_length:t,1),...
                'Ydata',Rn_ref(t-tail_length:t,2));
            set(tail_est,...
                'Xdata',Rn_est(t-tail_length:t,1),...
                'Ydata',Rn_est(t-tail_length:t,2));
        else
            set(tail_ref,...
                'Xdata',Rn_ref(1:t,1),...
                'Ydata',Rn_ref(1:t,2));
            set(tail_est,...
                'Xdata',Rn_est(1:t,1),...
                'Ydata',Rn_est(1:t,2));
        end        
        
        drawnow;
    end

end


%%
t = (1:Nsim) * dt;
figure('Name', 'Координаты');
subplot(2, 1, 1);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t(1:end-1), Rn_est(1:end-1, 1:2), '--', 'linewidth', 2);
plot(t(1:end-1), Rn_ref(1:end-1, 1:2), '-', 'linewidth', 1);
legend('X_{est}', 'Y_{est}', 'X_{ref}', 'Y_{ref}', 'fontsize', 16);
title('Координаты в навигационной СК, м', 'fontsize', 16);
ylabel('Координата, м', 'fontsize', 16);
xlabel('Время, сек', 'fontsize', 16);

subplot(2, 1, 2);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t, Rn_err(:, 1:2), 'linewidth', 1);
legend('\Delta X','\Delta Y', 'fontsize', 16);
title('Ошибки оценивания координат, м');
ylabel('Ошибка, м');
xlabel('Время, сек');

figure('Name','Ориентация');
subplot(2,1,1);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t(1:end-1), Euler_est(1:end-1, 1), '--', 'linewidth', 2);
plot(t(1:end-1), Euler_ref(1:end-1, 1),'-', 'linewidth', 1);
legend('\Theta_{est}', '\Theta_{ref}', 'fontsize', 16);
title('Угол курса', 'fontsize', 16);
ylabel('Угол, рад', 'fontsize', 16);
xlabel('Время, сек', 'fontsize', 16);

subplot(2,1,2);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t,Euler_err(:, 1), 'linewidth', 1);
legend('\Delta \Theta', 'fontsize', 16);
title('Ошибка оценивания угла курса', 'fontsize', 16);
ylabel('Ошибка, рад', 'fontsize', 16);
xlabel('Время, сек', 'fontsize', 16);

end





