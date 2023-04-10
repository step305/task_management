function fltr = imm_tracker(fltr, rn_platform, Cbn_platform, range, azimuth, update)
% Трекер координат, скоростей и модели движения подвижно-неподвижного
% препятствия
%
%   Входные аргументы:
%   fltr - структура IMM фильтра
%   rn_platform -  координаты платформы
%   Cbn_platform - углы ориентации платформы
%   range - дальность до препятствия
%   azimuth - угол пеленга препятствия
%   update - индикатор включения режима коррекции
%
%   Выходные аргументы:
%   fltr - структура IMM фильтра

% Параметры трекера
dt = 1.0;
q = [1e-2, 1e-4];
sigma_range = 0.2;
sigma_azimuth = 0.05;


%%
% Количество моделей движения
nmodels = 2;

%%
% Максимальная размерность вектора состояния
ndims = 4;

%%
% Индексы элементов вектора состояния
ind1 = [1 2 3 4];
ind2 = [1 2];

%%
% Цепь Маркова для вероятностей переключения моделей
prob = [
    0.9 0.1
    0.1 0.9
    ];

%%
% Предсказание
%
% Нормализующие коэффициенты для смешивания вероятностей моделей движения
c_j = zeros(1, nmodels);
for j = 1:nmodels
    for i = 1:nmodels
        c_j(j) = c_j(j) + prob(i, j).*fltr.mu(i);
    end
end

%%
% Смешивание вероятностей моделей движения
mu_ij = zeros(nmodels);
for i = 1:nmodels
    for j = 1:nmodels
        mu_ij(i, j) = prob(i, j) * fltr.mu(i) / c_j(j);
    end
end

%%
% Смешанные оценки вектора состояния для каждой модели движения
X_0j = zeros(ndims, nmodels);
for j = 1:nmodels
    X_0j(:, j) = zeros(ndims, 1);
    X_0j(ind1, j) = X_0j(ind1, j) + fltr.x1 * mu_ij(1, j);
    X_0j(ind2, j) = X_0j(ind2, j) + fltr.x2 * mu_ij(2, j);
end

%%
% Смешанные оценки матрицы ковариации для каждой модели движения
P_0j = zeros(nmodels,ndims,ndims);
for j = 1:nmodels
    P_0j(j,ind1,ind1) = squeeze(P_0j(j,ind1,ind1)) + mu_ij(1,j)*(fltr.P1 + (fltr.x1 - X_0j(ind1,j)) * (fltr.x1 - X_0j(ind1,j))');
    P_0j(j,ind2,ind2) = squeeze(P_0j(j,ind2,ind2)) + mu_ij(2,j)*(fltr.P2 + (fltr.x2 - X_0j(ind2,j)) * (fltr.x2 - X_0j(ind2,j))');
end

%%
% Предсказание
[fltr.x1, fltr.P1] = imm_predict(X_0j(ind1, 1), squeeze(P_0j(1, ind1, ind1)), 1, dt, q);
[fltr.x2, fltr.P2] = imm_predict(X_0j(ind2, 2), squeeze(P_0j(2, ind2, ind2)), 2, dt, q);


%%
% Коррекция
if update == 1
    lh = zeros(1, nmodels);

    %%
    % Корекция оценки для каждой модели движения
    [fltr.x1, fltr.P1, lh(1, 1)] = imm_update(fltr.x1, fltr.P1, ...
        rn_platform, Cbn_platform, ...
        range, azimuth, ...
        1, sigma_range, sigma_azimuth);

    [fltr.x2, fltr.P2, lh(1, 2)] = imm_update(fltr.x2, fltr.P2, ...
        rn_platform, Cbn_platform, ...
        range, azimuth, ...
        2, sigma_range, sigma_azimuth);

    %%
    % Вероятности моделей движения 
    c = sum(lh .* c_j);
    if (c > 0)
        fltr.mu = c_j .* lh/c;
    else
        fltr.mu = [0.5, 0.5];
    end
end

%%
% Обьединенные вектор состояния и матрица ковариации
X = zeros(ndims,1);
X(ind1, 1) = X(ind1, 1) + fltr.mu(1) * fltr.x1;
X(ind2, 1) = X(ind2, 1) + fltr.mu(2) * fltr.x2;

P = zeros(ndims);
P(ind1, ind1) = P(ind1, ind1) + fltr.mu(1) * (fltr.P1 + (fltr.x1 - X(ind1, 1)) * (fltr.x1 - X(ind1, 1))');
P(ind2, ind2) = P(ind2, ind2) + fltr.mu(2) * (fltr.P2 + (fltr.x2 - X(ind2, 1)) * (fltr.x2 - X(ind2, 1))');

fltr.x = X;
fltr.P = P;

end

%% Предсказание
function [x, P] = imm_predict(x, P, model, dt, q)
switch model
    case 1
        % Модель кинематики с дискретным белым шумом ускорений
        F = [
            1, 0, dt,   0
            0, 1,  0,  dt
            0, 0,  1,   0
            0, 0,  0,   1
            ];
        q1 = q(1);
        q2 = q(1);
        Q = [
            (dt^3*q1)/3,    0,          (dt^2*q1)/2, 0
            0,             (dt^3*q2)/3, 0,           (dt^2*q2)/2
            (dt^2*q1)/2,    0,          dt*q1,       0
            0,             (dt^2*q2)/2, 0,           dt*q2
            ];

        x = F * x;
        P = F * P * F' + Q;

    case 2
        % Модель кинематики неподвижного обьекта
        q1 = q(2);
        q2 = q(2);
        Q = [
            dt * q1, 0
            0, dt * q2
            ];
        P = P + Q;
end
end

%% Коррекция
function [x, P, lh] = imm_update(x, P, rn_platform, Cbn_platform, range, azimuth, model, sigma_range, sigma_azimuth)

Cnb_platform = Cbn_platform';

R = diag([sigma_range^2, sigma_azimuth^2]);

switch model
    case 1
        nstate = 4;
    case 2
        nstate = 2;
end

v = zeros(2, 1);
H = zeros(2, nstate);
I = eye(nstate);

xn = x;
xh = x;
lh = 0;
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

chi2 = (v' / S) * v;
if (chi2 < 1000)
    x = xn_;
    P = (I - K * H) * P * (I - K * H)' + K * R * K';
    lh = lh+studentLogprob([0, 0], S, 100.0, v');
end

lh = exp(lh);
end


function logp = studentLogprob(mu, Sigma, nu, x)
d = size(Sigma, 1);
x = bsxfun(@minus, x, rowvec(mu));
mahal = sum((x/Sigma).*x,2);
logc = gammaln(nu/2 + d/2) - gammaln(nu/2) - 0.5*logdet(Sigma) - (d/2)*log(nu) - (d/2)*log(pi);
logp = logc - (nu + d) / 2 * log1p(mahal / nu);
end

function y = logdet(A)
U = chol(A);
y = 2*sum(log(diag(U)));
end

function x = rowvec(x)
x = x(:)';
end
