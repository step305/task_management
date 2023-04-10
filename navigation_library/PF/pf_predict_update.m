function fltr = pf_predict_update(fltr, upsilon, dtheta, Anchors, Ranges, update, dt)
% Алгоритм фильтра частиц
%
%   Входные аргументы:
%   fltr - структура фильтра частиц
%   upsilon - управляющее воздействие робота - скорость продольного
%   перемещения
%   dtheta - управляющее воздействие робота - угловая скорость курсового
%   разворота
%   Ranges - измерения расстояний до базовых станций системы UWB
%   Anchors - координаты базовых станций системы UWB
%   update - индикатор доступности измерения от приемника UWB
%   dt - шаг интегрирования
%
%   Выходные аргументы:
%   fltr - структура фильтра частиц

for p = 1:fltr.N
    % Предсказание
    xp = fltr.p(:, p);
    xp = state_predict(xp, upsilon, dtheta, dt);
    fltr.p(:, p) = xp;

    % Коррекция
    if update == 1
        H  = dh(xp, Anchors);
        yp = h(xp, Anchors);
        P  = inv(fltr.invQ + H' * fltr.invR * H);
        P  = (P + P') / 2;
        mu = P * (fltr.invQ * xp + H' * fltr.invR * (Ranges - yp + H * xp));

        gauss.mean = mu;
        gauss.n = size(mu, 1);
        gauss.cov = P;
        gauss.const = 1 / sqrt((2 * pi)^size(mu, 1) * det(P));
        gauss.invCov = inv(P);
        [U, T] = schur(P);
        gauss.UsqrtT = U * (T .^ 0.5);

        % Выборка частицы
        fltr.p(:, p) = drawSamples(gauss, 1);
        fltr.w(:, p) = 1.0 / fltr.N;   

        % Вектор измерений
        y_pi = h(fltr.p(:, p), Anchors);

        % Плотность вероятности измерений
        likelihood = density(fltr.v_d, Ranges - y_pi);

        % Плотность вероятности системы
        prior = density(fltr.w_d, fltr.p(:,p) - xp);
        proposal = density(gauss, fltr.p(:,p));

        % Обновление веса частицы
        fltr.w(p) = fltr.w(p) * likelihood * prior / proposal;
    end

end

    %% Динамика вектора состояния 
    function x = state_predict(x, upsilon, dtheta, dt)
        xy = x(1:2,:);
        theta = x(3,:);

        upsilon_noise = 1e0;
        dtheta_noise  = 1e-1;

        upsilon = upsilon + randn * upsilon_noise;
        dtheta = dtheta + randn * dtheta_noise;

        xy(1) = xy(1) + upsilon * cos(theta) * dt;
        xy(2) = xy(2) + upsilon * sin(theta) * dt;
        theta = theta + dtheta * dt;

        x(1:2,:) = xy;
        x(3,:) = theta;
    end

    %% Матрица измерений
    function H = dh(x, Anchors)
        tag = [x(1:2, 1); 0.0];
        H = zeros(4,3);
        for i=1:4
            anchor = Anchors(i, :)';
            H(i, 1) = (tag(1) - anchor(1)) / norm(tag - anchor);
            H(i, 2) = (tag(2) - anchor(2)) / norm(tag - anchor);
        end
    end

    %% Вектор измерений
    function y = h(x, Anchors)
        tag = [x(1:2, 1); 0.0];
        y = zeros(4, 1);
        for i=1:4
            anchor = Anchors(i, :)';
            y(i, 1) = norm(tag - anchor);
        end
    end
end