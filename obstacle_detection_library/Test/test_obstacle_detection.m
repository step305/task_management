function [Rn_targets_hat, Vn_targets_hat] = test_obstacle_detection(Rn_ref, Euler_ref, Rn_targets, Theta_targets, ntargets)
% Тестирует работу детектора подвижных и неподвижных препятствий
%
% [Rn_ref, Euler_ref, sample_rate] = reference_trajectory_obstacle_detection(400);
% load('target_trajectories.mat')
% [Rn_targets_est, Vn_targets_est] = test_obstacle_detection(Rn_ref, Euler_ref, Rn_targets, Theta_targets, ntargets);
%
%   Входные аргументы:
%   Rn_ref -  координаты платформы
%   Euler_ref -  углы ориентации платформы
%   Rn_targets - координаты препятствий
%   Theta_targets - углы курса препятствий
%   ntargets - количество препятствий
%
%   Выходные аргументы:
%   Rn_targets_est -  оценки координат препятствий
%   Vn_targets_est -  оценки скоростей препятствий 


%% Модель радара
scan_angle = pi / 2;
range_max = 100;
range_min = 0.1;
measurement_rate = 100;

%% Инициализация анимации
close all;
fig = figure(...
    'NumberTitle', 'off',...
    'Clipping', 'on',...
    'ToolBar', 'figure', ...
    'Resize',  'on'...
    );
set (fig, 'Units', 'normalized', 'Position', [0,0,1,1]);
handle_axes_main = axes('position',[.25  .15  .5  .8]);
hold on, grid on;

set(handle_axes_main, 'Xlim', [-50 50]);
set(handle_axes_main, 'Ylim', [-50 50]);
set(handle_axes_main, 'Zlim', [-10 10]);
axis equal;
view(2);

% Оси опорного трехгранника платформы
at1 = [0; 0; 2.5]*5;
at2 = [0; 2.5; 0]*5;
at3 = [2.5; 0; 0]*5;
st =  [0; 0; 0];
pat1 = plot3([st(1) at1(1)],[st(2) at1(2)],...
    [st(3) at1(3)],'k--','linewidth', 2);
pat2 = plot3([st(1) at2(1)],[st(2) at2(2)],...
    [st(3) at2(3)],'k--','linewidth', 2);
pat3 = plot3([st(1) at3(1)],[st(2) at3(2)],...
    [st(3) at3(3)],'k--','linewidth', 2);

Vert = [ 0 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 0 0 1 ; 1 0 1 ; 1 1 1 ; 0 1 1];
Vert(:,1) = (Vert(:,1)-0.5)*7;
Vert(:,2) = (Vert(:,2)-0.5)*7;
Vert(:,3) = (Vert(:,3)-0.5)*7;
Faces = [ 1 2 6 5 ; 2 3 7 6 ; 3 4 8 7 ; 4 1 5 8 ; 1 2 3 4 ; 5 6 7 8 ];
ptch.Vertices = Vert;
ptch.Faces = Faces;
ptch.FaceVertexCData = jet(6);
ptch.FaceColor = 'flat';
patch_handle = patch(ptch);

% Опорные трехгранники подвижных препятствий
am1 = [0; 0; 5;];
am2 = [0; 5; 0;];
am3 = [5; 0; 0;];
pam1 = zeros(ntargets,1);
pam2 = zeros(ntargets,1);
pam3 = zeros(ntargets,1);
htargets = zeros(ntargets, 1);
cov_targets = zeros(ntargets, 1);
mu_targets = zeros(ntargets, 1);
for i=1:ntargets
    pam1(i,1) = plot3([0 am1(1)],[0 am1(2)], [0 am1(3)], '-','linewidth', 1, 'Color', 'k');
    pam2(i,1) = plot3([0 am2(1)],[0 am2(2)], [0 am2(3)], '-','linewidth', 1, 'Color', 'k');
    pam3(i,1) = plot3([0 am3(1)],[0 am3(2)], [0 am3(3)], '-','linewidth', 1, 'Color', 'k');
    htargets(i, 1) = plot(NaN, NaN, 'k.', 'MarkerSize', 20);
    cov_targets(i, 1) = plot(NaN,NaN, 'r', 'LineWidth', 1);
    mu_targets(i, 1) = plot(NaN,NaN, 'bo', 'LineWidth', 1);
end
handle_lines = plot(0, 0, '-.', 'Color', 'b', 'Linewidth', 0.5);
handle_measurements = plot(NaN, NaN, 'rs', 'MarkerSize', 7, 'LineWidth', 2);

% Траектории подвижных препятствий
for i=1:ntargets
    idx = 1:100:size(Rn_targets{i, 1}, 1);
    plot(Rn_targets{i, 1}(idx, 1), Rn_targets{i, 1}(idx, 2),':', 'Color', 'w', 'Linewidth', 0.1);
end

%% Структура детектора препятствий
% Структура трекера на базе EKF - для отслеживания подвижных препятствий
fltr_ekf.x  = [0; 0; 0; 0];
fltr_ekf.P  = diag([100, 100, 1, 1]);
fltr_ekf.hits = 0;
fltr_ekf.misses = 0;
fltr_ekf.latched = 0;
fltr_ekf.type = 'ekf';

% Структура трекера на базе IMM - для отслеживания подвижно-неподвижных
% препятствий
fltr_imm.x  = [0; 0; 0; 0];
fltr_imm.P  = diag([100, 100, 1, 1]);
fltr_imm.x1 = [0; 0; 0; 0];
fltr_imm.P1 = diag([100, 100, 1, 1]);
fltr_imm.x2 = [0; 0];
fltr_imm.P2 = diag([100, 100]);
fltr_imm.mu = [0.5, 0.5];
fltr_imm.hits = 0;
fltr_imm.misses = 0;
fltr_imm.latched = 0;  
fltr_imm.type = 'imm';

% Определение структуры детектора препятствий на базе разных трекеров
detector = {fltr_ekf, fltr_ekf, fltr_imm, fltr_imm, fltr_imm, fltr_imm, fltr_imm};

% Измерения радара - массивы
targets_range = zeros(ntargets, 1);
targets_azimuth = zeros(ntargets, 1);
targets_visible = zeros(ntargets, 1);
targets_pos = zeros(ntargets, 2);    

%% Цикл моделирования
Nsim = size(Rn_ref, 1);
Rn_targets_hat = {NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3)};
Vn_targets_hat = {NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3), NaN(Nsim, 3)};
for t=1:Nsim
    % Координаты платформы
    Rn_platform = Rn_ref(t, :)';
    Cbn_platform = angle_to_dcm(Euler_ref(t, 1), 0, 0)';

    % Измерения радара
    if (mod(t, measurement_rate) == 0)
        for j=1:ntargets
            track_length = size(Rn_targets{j, 1}, 1);
            if (t < track_length)
                Rn_target = Rn_targets{j, 1}(t, :)';
            else
                Rn_target = Rn_targets{j, 1}(track_length, :)';
            end
            [target_pos, range, azimuth, visible] = radar_model(Rn_platform, Cbn_platform, Rn_target, scan_angle, range_max, range_min);
            targets_pos(j, :) = target_pos;
            targets_range(j, 1) = range;
            targets_azimuth(j, 1) = azimuth;
            targets_visible(j, 1) = visible;
        end

        % Детектирование препятствий        
        detector = basic_obstacle_detector(detector, Rn_platform, Cbn_platform, targets_range, targets_azimuth, targets_visible);        
    end

    % Анимация
    if t < 1000
        animation_rate = 10;
    else
        animation_rate = 100;
    end

    if (mod(t, animation_rate) == 0)
        % Платформа
        a1_ = Cbn_platform*at1+Rn_platform;
        a2_ = Cbn_platform*at2+Rn_platform;
        a3_ = Cbn_platform*at3+Rn_platform;
        s_ = Rn_platform;
        set(pat1,'Xdata',[s_(1) a1_(1)]);
        set(pat1,'Ydata',[s_(2) a1_(2)]);
        set(pat1,'Zdata',[s_(3) a1_(3)]);
        set(pat2,'Xdata',[s_(1) a2_(1)]);
        set(pat2,'Ydata',[s_(2) a2_(2)]);
        set(pat2,'Zdata',[s_(3) a2_(3)]);
        set(pat3,'Xdata',[s_(1) a3_(1)]);
        set(pat3,'Ydata',[s_(2) a3_(2)]);
        set(pat3,'Zdata',[s_(3) a3_(3)]);
        Vert_ = Vert;
        for j=1:size(Vert,1)
            Vert_(j,:) = Vert(j,:) * Cbn_platform' + Rn_platform';
        end
        set(patch_handle,'Vertices',Vert_);

        % Подвижные препятствия и измерения
        z = [];
        k = 0;
        for i=1:ntargets
            track_length = size(Rn_targets{i,1}, 1);
            if (t < track_length)
                Rn_target = Rn_targets{i, 1}(t, :)';
                Cnb_target = angle_to_dcm(Theta_targets{i, 1}(t, 1), 0.0, 0.0);
            else
                Rn_target = Rn_targets{i, 1}(track_length, :)';
                Cnb_target = angle_to_dcm(Theta_targets{i, 1}(track_length, 1), 0.0, 0.0);
            end
            move_frame(Cnb_target, [Rn_target(1, 1); Rn_target(2, 1); Rn_target(3, 1)], pam1(i,1), pam2(i,1), pam3(i,1), am1, am2, am3);
            set(htargets(i, 1),  'Xdata', Rn_target(1, 1), 'Ydata', Rn_target(2, 1));

            % Измерения
            if (mod(t, measurement_rate) == 0)
                if (targets_visible(i, 1) == 1)
                    k = k + 1;
                    dX = (Rn_platform(1, 1) - targets_pos(i, 1));
                    dY = (Rn_platform(2, 1) - targets_pos(i, 2));
                    z(k,:) = [0 -dX -dY]; %#ok<AGROW>
                end
            end

            % Оценки координат препятствий
            if detector{i}.latched == 1
                p = covariance_ellipse(detector{i}.x(1:2, 1), detector{i}.P(1:2, 1:2));
                set(cov_targets(i, 1), 'xdata', p(1,:), 'ydata', p(2,:));
                set(mu_targets(i, 1), 'xdata', detector{i}.x(1, 1), 'ydata', detector{i}.x(2, 1));
                Rn_targets_hat{i}(t, :) = [t, detector{i}.x(1:2, 1)'];
                Vn_targets_hat{i}(t, :) = [t, detector{i}.x(3:4, 1)'];
            else
                set(cov_targets(i, 1), 'xdata', NaN, 'ydata', NaN);
                set(mu_targets(i, 1), 'xdata', NaN, 'ydata', NaN);
            end

        end

        set(handle_measurements, 'XData', targets_pos(targets_visible == 1, 1), 'YData', targets_pos(targets_visible == 1, 2));

        if ~isempty(z)
            plines = make_ranging_lines(z, [Rn_platform(1, 1), Rn_platform(2, 1)], size(z, 1));
            set(handle_lines, 'xdata', plines(1,:), 'ydata', plines(2,:));
        else
            set(handle_lines, 'xdata', 0, 'ydata', 0);
        end

        drawnow;
    end
end

for j=1:ntargets
    Rn_targets_hat{j} = rmmissing(Rn_targets_hat{j}(:, :));
    Vn_targets_hat{j} = rmmissing(Vn_targets_hat{j}(:, :));
end


end


