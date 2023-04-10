function p = covariance_ellipse(mu, Sigma)
% Вспомогательная функция для отображения эллипсоида ковариации

[U, D] = eig(Sigma);
n = 100;
t = linspace(0, 2*pi, n);
xy = [cos(t); sin(t)];
k = sqrt(6);
w = (k * U * sqrt(D)) * xy;
p = repmat(mu, [1 n]) + w;
end