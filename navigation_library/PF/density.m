function d = density(gDistr, x)
% Плотность вероятности гауссовкой случайной величины 

[nrow, ncol] = size(x);

if(nrow ~= gDistr.n)
  error('wrong input size for density calculation.');
end

d = zeros(1, ncol);
x = x - repmat(gDistr.mean, 1, ncol);


for i = 1:ncol
  d(i) = gDistr.const * exp(-x(:, i)' * gDistr.invCov * x(:, i) / 2);
end


  

