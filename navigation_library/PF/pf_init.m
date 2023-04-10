function fltr = pf_init(initPos, initAng)
% Инициализация фильтра частиц
%
%   Входные аргументы:
%   initPos - начальные коррдинаты робота в навигационной СК
%   initAng - начальный угол курса робота в навигационной СК
%
%   Выходные аргументы:
%   fltr - структура фильтра частиц

%%
% Количество частиц в фильтре
nParticle = 1000;

% Шумы системы
wPos = 1e-4;
wAng = 1e-3;

% Шумы измерений
wRnz = 1e0;

% Шум начальной выставки
wPosi = 1e2;
wAngi = 1e0;

% Шумы управлений
wUps = 1e-1;
wdThe = 1e-2;

% Размерности
nstate = 3;
nmeas = 4;

% Шумы системы
covW = diag([wPos, wPos, wAng]);
fltr.w_d.mean = zeros(nstate,1);
fltr.w_d.n = nstate;
fltr.w_d.cov = covW;
fltr.w_d.const = 1 / sqrt((2 * pi)^nstate * det(covW));
fltr.w_d.invCov = inv(covW);
[U, T] = schur(covW);
fltr.w_d.UsqrtT = U * (T .^ 0.5);

% Шумы измерений
covV = diag([wRnz wRnz wRnz wRnz]);
fltr.v_d.mean = zeros(nmeas,1);
fltr.v_d.n = nmeas;
fltr.v_d.cov = covV;
fltr.v_d.const = 1 / sqrt((2 * pi)^nmeas * det(covV));
fltr.v_d.invCov = inv(covV);
[U, T] = schur(covV);
fltr.v_d.UsqrtT = U * (T .^ 0.5);

% Начальное состояние
initMean = [initPos; initAng;];
initCov = diag([wPosi, wPosi, wAngi]);
initDistr.mean = initMean;
initDistr.n = nstate;
initDistr.cov = initCov;
initDistr.const = 1 / sqrt((2 * pi)^nstate * det(initCov));
initDistr.invCov = inv(initCov);
[U, T] = schur(initCov);
initDistr.UsqrtT = U * (T .^ 0.5);

% Выборка частиц начального состояния
fltr.N = nParticle;
fltr.p = drawSamples(initDistr, nParticle);
fltr.w = ones(1, nParticle) / nParticle;   
fltr.invQ = inv(covW);
fltr.invR = inv(covV);
fltr.Q = diag([wUps, wdThe]);

end
