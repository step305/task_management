function [Euler_, bw_] = test_AHRS(dWb, Fb, Mb, q_ref, q_init, fn, mn, dt)
%  [Euler_, bw_] = test_ahrs_dcm(dWb, Fb, Mb, q_ref, q_init, dt)
%  Моделирует функционирование алгоритма курсовертикали
%
%   Входные аргументы:
%   dWb -  Измерения датчиков угловой скорости 
%   Fb  - Измерения акселерометра
%   Mb  - Измерения магнитометра
%   q_ref  - Опорный кватернион ориентации
%   q_init - Начальное значение опорного кватерниона
%   fn - Вектор ускорения свободного падения 
%   mn - Вектор магнитного поля
%   dt - Шаг интегрирования
%
%   Выходные аргументы:
%   Euler_  - Оценки углов Эйлера
%   bw_     - Оценки сдвигов нулей датчиков угловой скорости

%%
close all; clc;

%% 
rng('shuffle');

%% Анимация
figure;
view(3);
axis equal;
hold on; grid on;
set(gca,'Xlim',[-2 2]);
set(gca,'Ylim',[-2 2]);
set(gca,'Zlim',[-2 2]);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
set(gcf,'renderer','opengl');

%% Оси опорного трехгранника
at1 = [0; 0; 2.5]*4;
at2 = [0; 2.5; 0]*4;
at3 = [2.5; 0; 0]*4;
st =  [0; 0; 0]*3;
pat1 = plot3([st(1) at1(1)],[st(2) at1(2)], [st(3) at1(3)],'b--','linewidth',0.1);
pat2 = plot3([st(1) at2(1)],[st(2) at2(2)], [st(3) at2(3)],'g--','linewidth',0.1);
pat3 = plot3([st(1) at3(1)],[st(2) at3(2)], [st(3) at3(3)],'r--','linewidth',0.1);

%% Оси оценки ориентации
ai1 = [0; 0; 1.5;];
ai2 = [0; 1.5; 0;];
ai3 = [1.5; 0; 0;];
sz =  [0; 0; 0;];
pai1 = plot3([sz(1) ai1(1)],[sz(2) ai1(2)], [sz(3) ai1(3)],'b-','linewidth',3);
pai2 = plot3([sz(1) ai2(1)],[sz(2) ai2(2)], [sz(3) ai2(3)],'g-','linewidth',3);
pai3 = plot3([sz(1) ai3(1)],[sz(2) ai3(2)], [sz(3) ai3(3)],'r-','linewidth',3);

%% Оси навигационной СК
an1 = [0; 0; 10;];
an2 = [0; 10; 0;];
an3 = [10; 0; 0;];
sz =  [0; 0; 0;];
plot3([sz(1) an1(1)],[sz(2) an1(2)], [sz(3) an1(3)],'b-','linewidth',0.1);
plot3([sz(1) an2(1)],[sz(2) an2(2)], [sz(3) an2(3)],'g-','linewidth',0.1);
plot3([sz(1) an3(1)],[sz(2) an3(2)], [sz(3) an3(3)],'r-','linewidth',0.1);
an1 = [0; 0; -10;];
an2 = [0; -10; 0;];
an3 = [-10; 0; 0;];
sz =  [0; 0; 0;];
plot3([sz(1) an1(1)],[sz(2) an1(2)], [sz(3) an1(3)],'b-','linewidth',0.1);
plot3([sz(1) an2(1)],[sz(2) an2(2)], [sz(3) an2(3)],'g-','linewidth',0.1);
plot3([sz(1) an3(1)],[sz(2) an3(2)], [sz(3) an3(3)],'r-','linewidth',0.1);

%%
Vert = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
Vert(:,1) = (Vert(:,1)-0.5);
Vert(:,2) = (Vert(:,2)-0.5);
Vert(:,3) = (Vert(:,3)-0.5);
Faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
ptch.Vertices = Vert;
ptch.Faces = Faces;
ptch.FaceVertexCData = copper(6);
ptch.FaceColor = 'flat';
patch_handle = patch(ptch);

%% 
th = text(-7,1,'Time ');
set(th,'FontSize',9,'Color','b','FontWeight','bold');

%% Начальные значения
% Случайные ошибки начальной ориентации
Cbn = quat_to_dcm(q_init);
err_x = randn*sqrt(2e-1);
err_y = randn*sqrt(2e-1);
err_z = randn*sqrt(2e-1);
Cbn = (eye(3)-skew([err_x; err_y; err_z]))*Cbn;

% Оценки сдвигов нулей ДУС
bw = [0; 0; 0];

% Матрица ковариации
P=zeros(6,6);
P(1:3,1:3)=diag([1e-2, 1e-2, 1e-2]); 
P(4:6,4:6)=diag([1e-6, 1e-6, 1e-6]); 

%% Логи
Nsim = size(dWb,1);  
Euler_ = zeros(Nsim,3);
bw_ = zeros(Nsim,3);
err_ = zeros(Nsim,3);

%% Цикл моделирования 
for i=1:Nsim
    
    
    %% Показания датчиков 
    dwb = dWb(i,:)';
    fb  = Fb(i,:)';
    mb  = Mb(i,:)';
    
    %% Алгоритм курсовертикали
    [Cbn, P, bw] = AHRS(Cbn, P, bw, dwb, fb, mb, fn, mn, dt);
    
    %% Сбор логов
    bw_(i,:) = bw*dt;
    [Euler_(i,1), Euler_(i,2), Euler_(i,3)] = dcm_to_angle(Cbn');
    
    % Ошибки углов ориентации
    Cbn_ref = quat_to_dcm(q_ref(i,:));
    Cerr = Cbn_ref*Cbn';
    [err_(i,1), err_(i,2), err_(i,3)] = dcm_to_angle(Cerr);
    


    
    %% Анимация
    if i < 300
        animation_rate = 1;
    else
        animation_rate = 30;
    end

    if(mod(i,animation_rate) == 0)
        
        %% Оси опорного трехгранника
        Cbnref = quat_to_dcm(q_ref(i,:))';
        a1_ = Cbnref'*at1;
        a2_ = Cbnref'*at2;
        a3_ = Cbnref'*at3;
        s_ = [0,0,0];
        set(pat1,'Xdata',[s_(1) a1_(1)]);
        set(pat1,'Ydata',[s_(2) a1_(2)]);
        set(pat1,'Zdata',[s_(3) a1_(3)]);
        set(pat2,'Xdata',[s_(1) a2_(1)]);
        set(pat2,'Ydata',[s_(2) a2_(2)]);
        set(pat2,'Zdata',[s_(3) a2_(3)]);
        set(pat3,'Xdata',[s_(1) a3_(1)]);
        set(pat3,'Ydata',[s_(2) a3_(2)]);
        set(pat3,'Zdata',[s_(3) a3_(3)]);
        
        %% Оценка ориентации трехгранника 
        a1_ = Cbn*ai1;
        a2_ = Cbn*ai2;
        a3_ = Cbn*ai3;
        s_ = [0, 0, 0];
        set(pai1,'Xdata',[s_(1) a1_(1)]);
        set(pai1,'Ydata',[s_(2) a1_(2)]);
        set(pai1,'Zdata',[s_(3) a1_(3)]);
        set(pai2,'Xdata',[s_(1) a2_(1)]);
        set(pai2,'Ydata',[s_(2) a2_(2)]);
        set(pai2,'Zdata',[s_(3) a2_(3)]);
        set(pai3,'Xdata',[s_(1) a3_(1)]);
        set(pai3,'Ydata',[s_(2) a3_(2)]);
        set(pai3,'Zdata',[s_(3) a3_(3)]);
        
        %
        Vert_ = Vert;
        for j=1:size(Vert,1)
            Vert_(j,:) = (Vert(j,:)*Cbn');
        end
        set(patch_handle,'Vertices',Vert_);
        
        set(th,'String',...
            sprintf('tm %2.1f\nps   % 7.5f\ntt     % 7.5f\ngm  % 7.5f',...
            i*dt,Euler_(i,1)*180/pi,Euler_(i,2)*180/pi,Euler_(i,3)*180/pi));
        
        drawnow;
    end
end

%% 
t = (1:Nsim)*dt;

figure('Name','Ошибки углов ориентации');
hold on;
grid on;
plot(t,err_*180/pi);
legend('\delta\psi','\delta\theta','\delta\gamma');
ylabel('Angles errors, deg');
xlabel('Time,sec');

figure('Name','Оценки сдвигов нулей ДУС');
hold on;
grid on;
plot(t,bw_);
legend('\delta\omega_x','\delta\omega_y','\delta\omega_z');
ylabel('Gyro Bias, rad/sec');
xlabel('Time,sec');


end

