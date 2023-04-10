function X = drawSamples(gDistr, totalNum)
% Выборка из гауссовского распределения 

X = gDistr.UsqrtT * randn(gDistr.n, totalNum) + repmat(gDistr.mean, 1, totalNum);
