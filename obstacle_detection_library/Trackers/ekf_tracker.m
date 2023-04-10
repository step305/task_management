function fltr = ekf_tracker(fltr, rn_platform, Cbn_platform, range, azimuth, update)
% Трекер координат и скоростей подвижного препятствия
%
%   Входные аргументы:
%   fltr - структура EKF фильтра
%   rn_platform -  координаты платформы
%   Cbn_platform - углы ориентации платформы
%   range - дальность до препятствия
%   azimuth - угол пеленга препятствия
%   update - индикатор включения режима коррекции
%
%   Выходные аргументы:
%   fltr - структура EKF фильтра


% Параметры трекера
dt = 1.0;
q = 1e-2;
sigma_range = 0.2;
sigma_azimuth = 0.05;


x = fltr.x;
P = fltr.P;

% Предсказание
% Модель кинематики с дискретным белым шумом ускорений
F = [
    1, 0, dt, 0
    0, 1, 0,  dt
    0, 0, 1,  0
    0, 0, 0,  1
    ];

Q = [
    (dt^3*q)/3,    0,          (dt^2*q)/2, 0
    0,             (dt^3*q)/3, 0,           (dt^2*q)/2
    (dt^2*q)/2,    0,          dt*q,       0
    0,             (dt^2*q)/2, 0,           dt*q
    ];

x = F * x ;
P = F * P * F' + Q;

% Коррекция
if update == 1
    Cnb_platform = Cbn_platform';
    R = diag([sigma_range^2, sigma_azimuth^2]);

    v = zeros(2, 1);
    H = zeros(2, 4);
    I = eye(4);

    xn = x;
    xh = x;
    niter = 10;
    for n=1:niter
        rn_target_hat = [xn(1:2, 1); 0];

        fb = Cnb_platform * (rn_target_hat - rn_platform);
        azimuth_hat = atan2(fb(2, 1), fb(1, 1));

        range_hat = norm(rn_target_hat - rn_platform);

        v(1, 1) = range - range_hat;
        v(2, 1) = pi2pi(azimuth - azimuth_hat);

        H(1, 1) = (rn_target_hat(1, 1) - rn_platform(1, 1)) / norm(rn_target_hat - rn_platform);
        H(1, 2) = (rn_target_hat(2, 1) - rn_platform(2, 1)) / norm(rn_target_hat - rn_platform);

        c1_1 = Cnb_platform(1, 1);
        c1_2 = Cnb_platform(1, 2);
        c1_3 = Cnb_platform(1, 3);
        c2_1 = Cnb_platform(2, 1);
        c2_2 = Cnb_platform(2, 2);
        c2_3 = Cnb_platform(2, 3);

        rp1 = rn_platform(1, 1);
        rp2 = rn_platform(2, 1);
        rp3 = rn_platform(3, 1);

        rt1 = rn_target_hat(1, 1);
        rt2 = rn_target_hat(2, 1);
        rt3 = rn_target_hat(3, 1);
        
        H(2, 1) = -(c2_1*(c1_2*(rp2 - rt2) + c1_3*(rp3 - rt3)) - c1_1*(c2_2*(rp2 - rt2) + ...
            c2_3*(rp3 - rt3)))/((c1_1*rp1 + c1_2*rp2 + c1_3*rp3 - c1_1*rt1 - c1_2*rt2 - ...
            c1_3*rt3)^2 + (c2_1*rp1 + c2_2*rp2 + c2_3*rp3 - c2_1*rt1 - c2_2*rt2 - c2_3*rt3)^2);

        H(2, 2) = -(c2_2*(c1_1*(rp1 - rt1) + c1_3*(rp3 - rt3)) - c1_2*(c2_1*(rp1 - rt1) + ...
            c2_3*(rp3 - rt3)))/((c1_1*rp1 + c1_2*rp2 + c1_3*rp3 - c1_1*rt1 - c1_2*rt2 - ...
            c1_3*rt3)^2 + (c2_1*rp1 + c2_2*rp2 + c2_3*rp3 - c2_1*rt1 - c2_2*rt2 - c2_3*rt3)^2);

        S = H * P * H' + R;
        K = (P * H') / S;
        xn_ = xh + K * (v - H * (xh - xn));

        if norm(xn_ - xn) < 1e-7
            break;
        end
        xn = xn_;
    end

    x = xn_;
    P = (I - K * H) * P * (I - K * H)' + K * R * K';

end

fltr.x = x;
fltr.P = P;

end


