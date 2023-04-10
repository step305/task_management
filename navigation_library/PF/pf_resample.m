function fltr = pf_resample(fltr, threshold)
% Ресэмплинг частиц

Nmin = fltr.N * threshold;
sum_w = sum(fltr.w);
fltr.w = fltr.w / sum_w;

[keep, Neff] = stratified_resample(fltr.w);
if Neff < Nmin
    fltr.p(:, 1:fltr.N) = fltr.p(:,keep);
    fltr.w(:, :) = ones(1, fltr.N) / fltr.N;
end

%%
function [keep, Neff] = stratified_resample(w)

w = w / sum(w); 
Neff = 1 / sum(w .^ 2);

len = length(w);
keep = zeros(1,len);
select = stratified_random(len);
w = cumsum(w);

ctr = 1;
for i = 1:len
    while ctr<=len && select(ctr)<w(i)
        keep(ctr) = i;
        ctr = ctr+1;
    end
end
end

%% 
function s = stratified_random(N)

k = 1/N;
di = (k/2):k:(1-k/2); 
s = di + rand(1,N) * k - (k/2);
end

end