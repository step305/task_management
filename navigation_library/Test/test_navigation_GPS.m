function [Rn_, Vb_, Euler_, Rn_err, Vb_err, Euler_err] = test_navigation_GPS(Rn_ref, Vb_ref, Euler_ref, dThe_ref, dt)
%  Моделирует функционирование алгоритма системы навигации колесного робота
%  с коррекцией от GPS - измерения GPS предварительно спроецированны в
%  локальную горизонтальную систему координат
%
% [dThe_ref, dUps_ref, Rn_ref, Vb_ref, Euler_ref, Lat_ref, Lon_ref, dThet_ref, sample_rate] = reference_trajectory(600);
% [Rn_, Vb_, Euler_, Rn_err, Vb_err, Euler_err] = test_navigation_GPS(Rn_ref, Vb_ref, Euler_ref, dThe_ref, sample_rate);
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
%   Vb_ - Оценки значений скоростей в связанной СК
%   Euler_ - Оценки значений углов ориентации
%   Rn_err - Ошибки оценок значений координат в навигационной СК
%   Vb_err - Ошибки оценок значений скоростей в связанной СК
%   Euler_err - Ошибки оценок значений углов ориентации

close all;

%% 
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
ptch.FaceVertexCData = copper(6);
ptch.FaceColor = 'flat';
patch_handle = patch(ptch);

%%
hu = plot3(0,0,-3,'k.','MarkerSize', 20);
ht = plot(0,0,'ks','MarkerSize',10);
th = text(-55,-55, 'Time ');
set(th,'FontSize',12,'Color','b');
tail_ref = plot(0,0,'k.-');
tail_est = plot(0,0,'b.-');
tail_length = 1800;

%% Параметры GPS
gps_pos_noise = 5e0;
gps_vel_noise = 1e-1;
gps_lever = [1.0; 1.0; 0.0]; 
gps_rate  = 100;

%% Параметры одометра
lr = 0.3; 
ll = 0.3; 
d  = 0.4; 
err_lr = -2e-2;
err_ll =  2e-2;
odometer_noise = 1e-5;

%% Gyro Parameters
gyro_bias = 1e-3;
gyro_noise = 1e-5;

%% Начальные условия
rn = [Rn_ref(1,1) + 20.0;
    Rn_ref(1,2) - 20.0;
    0];
Cbn = angle_to_dcm(deg2rad(90), 0, 0)';

%% Мартрица ковариации 
P=zeros(6);
P(1:2,1:2)=eye(2) * 1e0;  %position errors
P(3,3)=1e0;               %heading error
P(4:5,4:5)=eye(2) * 1e-3; %wheel radiuses errors
P(6,6)=1e-5;              %gyro bias  error

% Оценки ошибок датчиков
dlr = 0;
dll = 0;
dw  = 0;

%% Логи
Nsim = size(Rn_ref, 1)-1; 
Rn_ = zeros(Nsim,3);
Vb_ = zeros(Nsim,3);
Euler_ = zeros(Nsim,3);
Euler_err = zeros(Nsim,3);
Rn_err = zeros(Nsim,3);
Vb_err = zeros(Nsim,3);
Vx_odo = zeros(Nsim,1);
Wz_odo = zeros(Nsim,1);
dlr_ = zeros(Nsim,1);
dll_ = zeros(Nsim,1);
dw_  = zeros(Nsim,1);
Pos_ref = zeros(Nsim,2);
Pos_est = zeros(Nsim,2);
cnt = 0;

%% Цикл моделирования
for i=2:Nsim

    %% Опорные измерения
    Rnref = Rn_ref(i, :)';
    Cnb_ref = angle_to_dcm(Euler_ref(i, 1), Euler_ref(i, 2), Euler_ref(i, 3));
    Cbn_ref = Cnb_ref';
    Vbref = Vb_ref(i, :)';
    Wbref = dThe_ref(i, :)'/dt;

    %% Измерения GPS 
    rg_ref = Rnref + Cbn_ref * gps_lever;
    pos_gps = rg_ref + gps_pos_noise * randn(3, 1);
    vg_ref = Cbn_ref * (Vbref + skew(Wbref) * gps_lever);
    vel_gps = vg_ref + gps_vel_noise * randn(3, 1);
    if (mod(i, gps_rate) == 0)
        gps_update = true;
    else
        gps_update = false;
    end 

    %% Одометр
    fir = Vb_ref(i, 1) / (lr + err_lr) + (dThe_ref(i, 3) / dt) * d / (2 * (lr + err_lr)) + odometer_noise * randn;
    fil = Vb_ref(i, 1) / (ll + err_ll) - (dThe_ref(i, 3) / dt) * d / (2 * (ll + err_ll)) + odometer_noise * randn;

    %% Курсовой ДУС
    dThe = dThe_ref(i,3) + gyro_bias + gyro_noise * randn;

    %% Алгоритм навигации
    [rn, Cbn, Vb_odometer, Wb_odometer, P, dw, dlr, dll] = ...
        navigation_system_GPS(rn, Cbn, P, dt, dThe, fir, fil, lr, ll, d, ...
        pos_gps, vel_gps, gps_lever, gps_update, dw, dlr, dll);

    %% 
    Rn_(i-1,:) = rn;
    Vb_(i-1,:) = [Vb_odometer, 0, 0];
    [Euler_(i-1,1),~,~] = dcm_to_angle(Cbn');
    Vx_odo(i-1,1) = Vb_odometer;
    Wz_odo(i-1,1) = Wb_odometer;
    dlr_(i-1,:)   = dll;
    dll_(i-1,:)   = dlr;
    dw_(i-1,:)    = dw;
    
    % Ошибки
    Cref = angle_to_dcm(Euler_ref(i,1), Euler_ref(i,2), Euler_ref(i,3));
    Cerr = Cref*Cbn;
    [Euler_err(i-1,1), Euler_err(i-1,2), Euler_err(i-1,3)] = dcm_to_angle(Cerr);
    Rn_err(i-1,:) = Rn_ref(i,:) - rn';
    Vb_err(i-1,:) = Vb_ref(i,:)-[Vb_odometer, 0, 0];

    %% Анимация
    if i < 3000
        animation_rate = 10;
    else
        animation_rate = 100;
    end    
    if(mod(i,animation_rate) == 0)
        cnt = cnt+1;

        set(ht, 'Xdata',Rn_ref(i,1), 'Ydata',Rn_ref(i,2));
        set(hu, 'Xdata', pos_gps(1,1), 'Ydata', pos_gps(2,1));
        Pos_ref(cnt,:) = Rn_ref(i,1:2);

        Cbnref = angle_to_dcm(Euler_ref(i,1),Euler_ref(i,2),Euler_ref(i,3));
        a1_ = Cbnref'*at1+Rn_ref(i,:)';
        a2_ = Cbnref'*at2+Rn_ref(i,:)';
        a3_ = Cbnref'*at3+Rn_ref(i,:)';
        s_ = Rn_ref(i,:);
        set(pat1,'Xdata',[s_(1) a1_(1)]);
        set(pat1,'Ydata',[s_(2) a1_(2)]);
        set(pat1,'Zdata',[s_(3) a1_(3)]);
        set(pat2,'Xdata',[s_(1) a2_(1)]);
        set(pat2,'Ydata',[s_(2) a2_(2)]);
        set(pat2,'Zdata',[s_(3) a2_(3)]);
        set(pat3,'Xdata',[s_(1) a3_(1)]);
        set(pat3,'Ydata',[s_(2) a3_(2)]);
        set(pat3,'Zdata',[s_(3) a3_(3)]);

        Pos_est(cnt,:) = Rn_(i-1,1:2);
        a1_ = Cbn*ai1+Rn_(i-1,:)';
        a2_ = Cbn*ai2+Rn_(i-1,:)';
        a3_ = Cbn*ai3+Rn_(i-1,:)';
        s_ = Rn_(i-1,:);
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
            Vert_(j,:) = (Vert(j,:)*Cbn')+[Rn_(i-1,1) Rn_(i-1,2) Rn_(i-1,3)];
        end
        set(patch_handle,'Vertices',Vert_);

        set(th,'String',sprintf('Время %2.1f сек.',i*dt));

        if (cnt > tail_length)
            set(tail_ref,...
                'Xdata',Pos_ref(cnt-tail_length:cnt,1),...
                'Ydata',Pos_ref(cnt-tail_length:cnt,2));
            set(tail_est,...
                'Xdata',Pos_est(cnt-tail_length:cnt,1),...
                'Ydata',Pos_est(cnt-tail_length:cnt,2));
        else
            set(tail_ref,...
                'Xdata',Pos_ref(1:cnt,1),...
                'Ydata',Pos_ref(1:cnt,2));
            set(tail_est,...
                'Xdata',Pos_est(1:cnt,1),...
                'Ydata',Pos_est(1:cnt,2));
        end
        drawnow;
    end
end

%%
t = (1:Nsim)*dt;
figure('Name', 'Координаты');
subplot(2,1,1);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t(1:end-1), Rn_(1:end-1, 1:2), '--', 'linewidth', 2);
plot(t(1:end-1), Rn_ref(2:end-1, 1:2), '-', 'linewidth', 1);
legend('X_{est}', 'Y_{est}', 'X_{ref}', 'Y_{ref}', 'fontsize', 16);
title('Координаты в навигационной СК, м', 'fontsize', 16);
ylabel('Координата, м', 'fontsize', 16);
xlabel('Время, сек', 'fontsize', 16);

subplot(2,1,2);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t,Rn_err(:, 1:2), 'linewidth', 1);
legend('\Delta X','\Delta Y', 'fontsize', 16);
title('Ошибки оценивания координат, м');
ylabel('Ошибка, м');
xlabel('Время, сек');

figure('Name','Скорость');
subplot(2,1,1);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t(1:end-1),Vb_(1:end-1, 1), '--', 'linewidth', 2);
plot(t(1:end-1),Vb_ref(2:end-1, 1),'-', 'linewidth', 1);
legend('Vx_{est}', 'Vx_{ref}', 'fontsize', 16);
title('Скорости в связанной СК, м/c', 'fontsize', 16);
ylabel('Скорость, м/с', 'fontsize', 16);
xlabel('Время, сек', 'fontsize', 16);

subplot(2,1,2);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t,Vb_err(:, 1), 'linewidth', 1);
legend('\Delta Vx', 'fontsize', 16);
title('Ошибка оценивания скорости, м/c', 'fontsize', 16);
ylabel('Ошибка, м/с', 'fontsize', 16);
xlabel('Время, сек', 'fontsize', 16);

figure('Name','Ориентация');
subplot(2,1,1);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t(1:end-1), Euler_(1:end-1, 1), '--', 'linewidth', 2);
plot(t(1:end-1), Euler_ref(2:end-1, 1),'-', 'linewidth', 1);
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

figure('Name','Ошибки одометра');
subplot(2,1,1);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t(1:end-1), Vx_odo(1:end-1,1)-Vb_ref(2:end-1,1), 'linewidth', 1);
legend('\Delta Vx_{odo}', 'fontsize', 16);
title('Ошибки оценивания скорости одометром, м/c', 'fontsize', 16);
ylabel('Ошибка, м/с', 'fontsize', 16);
xlabel('Время, сек', 'fontsize', 16);

subplot(2,1,2);
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t, dll_, '--', 'linewidth', 2);
plot(t, ones(length(t),1)*err_ll, '-', 'linewidth', 1);
plot(t, dlr_, '--', 'linewidth', 2);
plot(t, ones(length(t),1)*err_lr, '-', 'linewidth', 1);
legend('\delta ll_{est}', '\delta ll_{ref}', '\delta lr_{est}', '\delta lr_{ref}', 'fontsize', 16);
title('Оценки погрешностей длин окружностей колес', 'fontsize', 16);
ylabel('Погрешность, м', 'fontsize', 16);
xlabel('Время, сек', 'fontsize', 16);

figure('Name','Ошибки ДУС');
ax = gca;
ax.FontSize = 14;
hold on; grid on;
plot(t(1:end-1), dw_(1:end-1, :)*dt, '--', 'linewidth', 2);
plot(t(1:end-1), ones(length(t)-1, 1)*gyro_bias, '-', 'linewidth', 1);
legend('\Delta \omega_{est}', '\Delta \omega_{ref}', 'fontsize', 16);
title('Оценка сдвига нуля ДУС');
ylabel('Сдвиг нуля, рад/с');
xlabel('Время, сек');



