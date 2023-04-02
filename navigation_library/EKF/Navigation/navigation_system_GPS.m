function [rn, Cbn, Vb_odometer, Wb_odometer, P, dw, dlr, dll] = ...
    navigation_system_GPS(rn, Cbn, P, dt, dThe, fir, fil, lr, ll, d, ...
    pos_gps, vel_gps, lb, gps_update, dw, dlr, dll)


%% Одометр
% компенсация ошибок одометра - длины правого и левого колес и длины оси
% колесной пары
lr = lr - dlr;
ll = ll - dll;

% Скорости продольного перемещения и курсового разворота
Vb_odometer = 1/2 * (fir * lr + fil * ll);
Wb_odometer = 1/d * (fir * lr - fil * ll);

% Компенсация сдвига нуля гироскопа
dthet_z = dThe - dw * dt;

% Обновление матрицы направляющих косинусов (интегрирование ориентации)
rot_norm = norm(dthet_z);
sr_a = 1 - (rot_norm ^ 2 / 6) + (rot_norm ^ 4 / 120);
sr_b = (1 / 2) - (rot_norm ^ 2 / 24) + (rot_norm ^ 4 / 720);
Cbb = [
    1 - sr_b * dthet_z ^ 2, -sr_a * dthet_z, 0;
    sr_a * dthet_z, 1 - sr_b * dthet_z ^ 2, 0;
    0, 0, 1
    ];
Cbn=Cbn * Cbb;

% Обновление координат
rn = rn + Cbn * [Vb_odometer; 0; 0] * dt;

%% Фильтр Калмана - предсказание

% Матрица динамики системы в непрерывном времени
An = zeros(3);
An(1,3) =  Cbn(2,1) * Vb_odometer;
An(2,3) = -Cbn(1,1) * Vb_odometer;

% Матрица динамики системы в дискретном времени
Fnd=eye(3) + An * dt + An * An * dt * dt / 2;

% Матрица влияния шумов системы
Nn = zeros(3, 3);
Nn(1, 1) =  Cbn(1, 1) * fir / 2;
Nn(1, 2) =  Cbn(1, 1) * fil / 2;
Nn(2, 1) =  Cbn(2, 1) * fir / 2;
Nn(2, 2) =  Cbn(2, 1) * fil / 2;
Nn(3, 3) = -1;

% Матрица ковариации шумов системы
pos_noise = 1e-6; 
att_noise = 1e-8; 
R = diag([pos_noise, pos_noise, att_noise]);
Qn = Nn * R * Nn';

% Дискретизация с интегрированием методом трапеций
Qnd = dt / 2 * (Fnd * Qn + Qn * Fnd');

% Матрица динамики системы и матрица ковариации шумов измерений 
F = zeros(6);
Q = zeros(6);
F(1:3, 1:3) = Fnd;
F(4:6, 4:6) = eye(3);
F(1:3, 4:6) = Nn * dt;

Q(1:3, 1:3) = Qnd;
nlr = 1e-12;
nll = 1e-12;
nw  = 1e-20;
Qwhe = diag([nlr, nll, nw]);
Q(4:6, 4:6) = Qwhe;
Q(1:3, 4:6) = Nn * Qwhe * dt / 2;
Q(4:6, 1:3) = Q(1:3, 4:6)';

% Обновление матрицы ковариации фильтра Калмана
P = EKF_prediction(P, F, Q);

%% Фильтр Калмана - коррекция, GPS
if (gps_update)
        % Вектор измерений
        v = zeros(4, 1);
        pos_gps_hat = rn + Cbn * lb;
        v(1:2, 1) = pos_gps_hat(1:2, 1) - pos_gps(1:2, 1);

        vb = [Vb_odometer; 0; 0];
        wb = [0; 0; Wb_odometer];
        vel_gps_hat = Cbn * (vb + skew(wb) * lb);
        v(3:4, 1) = vel_gps_hat(1:2, 1) - vel_gps(1:2, 1);   

        % Матрица измерений
        H = zeros(4, 6);
        H(1:2, 1:2) = eye(2);
        dpsi = skew(Cbn * lb);
        H(1:2, 3) = dpsi(1:2, 3);    
        dpsi = skew(vel_gps_hat);
        H(3:4, 3) =  dpsi(1:2, 3);
        domega   = -Cbn * skew(lb);
        H(3:4, 6) = domega(1:2, 3);
        H(3, 4)   = Cbn(1, 1) * fir / 2;
        H(3, 5)   = Cbn(1, 1) * fil / 2;
        H(4, 4)   = Cbn(2, 1) * fir / 2;
        H(4, 5)   = Cbn(2, 1) * fil / 2;

        % Шумы измерений
        R(1:2, 1:2) = eye(2) * 1e-1;
        R(3:4, 3:4) = eye(2) * 1e-2;

        % Обновление
        [x, P] = EKF_correction(P, v, H, R);

        % Коррекция оценки координат
        rn = rn - [x(1, 1); x(2, 1); 0];

        % Коррекция оценки ориентации
        En  = [1, -x(3, 1), 0; x(3, 1), 1, 0; 0, 0, 1];
        Cbn = En * Cbn;
        Cbn = dcmnormalize(Cbn);

        % Коррекция оценок ошибок датчиков
        dlr = dlr + x(4, 1);
        dll = dll + x(5, 1);
        dw  = dw + x(6, 1);
end

%%  Фильтр Калмана - коррекция, одометр
% Вектор измерений 
v = dthet_z / dt - Wb_odometer;

% Матрица измерений
H = zeros(1, 6);
H(1, 4) = -fir / d;
H(1, 5) =  fil / d;
H(1, 6) =  1;

% Шумы измерений
R = 1e-2;

% Обновление 
[x, P] = EKF_correction(P, v, H, R);

% Коррекция оценки координат
rn = rn - [x(1, 1); x(2, 1); 0];

% Коррекция оценки ориентации
En  = [1, -x(3, 1), 0; x(3, 1), 1, 0; 0, 0, 1];
Cbn = En * Cbn;
Cbn = dcmnormalize(Cbn);

% Коррекция оценок ошибок датчиков
dlr = dlr + x(4, 1);
dll = dll + x(5, 1);
dw  = dw + x(6, 1);

end