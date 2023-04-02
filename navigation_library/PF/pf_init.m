function fltr_ = pf_init(initPos, initAng)


%%
% Количество частиц фильтра
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
fltr_.w_d.mean = zeros(nstate,1);
fltr_.w_d.n = nstate;
fltr_.w_d.cov = covW;
fltr_.w_d.const = 1 / sqrt((2 * pi)^nstate * det(covW));
fltr_.w_d.invCov = inv(covW);
[U, T] = schur(covW);
fltr_.w_d.UsqrtT = U * (T .^ 0.5);

% Шумы измерений
covV = diag([wRnz wRnz wRnz wRnz]);
fltr_.v_d.mean = zeros(nmeas,1);
fltr_.v_d.n = nmeas;
fltr_.v_d.cov = covV;
fltr_.v_d.const = 1 / sqrt((2 * pi)^nmeas * det(covV));
fltr_.v_d.invCov = inv(covV);
[U, T] = schur(covV);
fltr_.v_d.UsqrtT = U * (T .^ 0.5);

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
fltr_.N = nParticle;
fltr_.p = drawSamples(initDistr, nParticle);
fltr_.w = ones(1, nParticle) / nParticle;   
fltr_.invQ = inv(covW);
fltr_.invR = inv(covV);
fltr_.Q = diag([wUps, wdThe]);

end
