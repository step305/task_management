function [rn, Cbn, Vb_odometer, Wb_odometer, P, dlu, dw, dlr, dll, dd] = ...
    navigation_system_UWB(rn, Cbn, P, dt, dThe, fir, fil, lr, ll, d, ...
    Ranges, Anchors, lu, uwb_update, dlu, dw, dlr, dll, dd)


%% UWB
% компенсация ошибок установки приемника UWB
lu  = lu - dlu;

%% Одометр
% компенсация ошибок одометра - длины правого и левого колес и длины оси
% колесной пары
lr = lr - dlr;
ll = ll - dll;
d  = d - dd;

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
Cnb=Cbn';

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
nlr = 1e-6; 
nll = 1e-6; 
nw  = 1e-8; 
R = diag([nlr, nll, nw]);
Qn = Nn * R * Nn';

% Дискретизация с интегрированием методом трапеций
Qnd = dt / 2 * (Fnd * Qn + Qn * Fnd');

% Матрица динамики системы и матрица ковариации шумов измерений 
F = zeros(9);
Q = zeros(9);
F(1:3, 1:3) = Fnd;
F(4:9, 4:9) = eye(6);
F(1:3, 4:6) = Nn * dt;

Q(1:3, 1:3) = Qnd;
nlr = 1e-12;
nll = 1e-12;
nw  = 1e-20;
Qwhe = diag([nlr, nll, nw]);
Q(4:6, 4:6) = Qwhe;
Q(1:3, 4:6) = Nn * Qwhe * dt / 2;
Q(4:6, 1:3) = Q(1:3, 4:6)';
nd  = 1e-12; 
Q(7, 7) = nd;
nlu = 1e-16; 
Qlev = diag([nlu, nlu]);
Q(8:9, 8:9) = Qlev;

% Обновление матрицы ковариации фильтра Калмана
P = EKF_prediction(P, F, Q);

%% Фильтр Калмана - коррекция, UWB
if (uwb_update)
    for j=1:size(Ranges, 1)
        % Вектор измерений
        l = [lu; 0];
        ru = rn + Cbn * [lu; 0];
        ra = Anchors(j, :)';
        r = norm(ru - ra);
        v = r - Ranges(j);

        % Матрица измерений
        H = zeros(1, 9);
        dr0 = rn' - ra' + l' * Cnb;
        H(1, 1:2) = 1 / r * dr0(1, 1:2);
        dl = l' + (rn' - ra') * Cbn;
        H(1,8:9) =  1 / r * dl(1, 1:2);
        dpsi = (rn' - ra' + l' * Cnb) * skew(Cbn * l);
        H(1,3) =  1 / r * dpsi(1, 3);

        % Шумы измерений
        R = 5e-2;

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
        dd  = dd + x(7, 1);
        dlu = dlu + x(8:9, 1);
    end
end

%%  Фильтр Калмана - коррекция, одометр
% Вектор измерений 
v = dthet_z / dt - Wb_odometer;

% Матрица измерений
H = zeros(1, 9);
H(1, 4) = -fir / d;
H(1, 5) =  fil / d;
H(1, 6) =  1;
H(1, 7) = (fir * lr - fil * ll) / d ^ 2;

% Шумы измерений
R = 1e-2;

% Обновление 
I = eye(9);
K = (P * H') / (H * P * H' + R);
P = (I - K * H) * P * (I - K * H)' + K * R * K';
x = K * v;

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
dd  = dd + x(7, 1);
dlu = dlu + x(8:9, 1);

end
