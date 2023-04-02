function X = drawSamples(gDistr, totalNum)
X = gDistr.UsqrtT * randn(gDistr.n, totalNum) + repmat(gDistr.mean, 1, totalNum);
