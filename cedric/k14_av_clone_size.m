function [b, s, t] = k14_av_clone_size(r, gamma, mu, theta, ts)

% For k14, we make the following assumption: theta is the proportion of
% induced clones which are stem like, and they immediatedly (possible with
% some time delay) divide to give another stem cell and a CP cell: S ->
% S+A, however new stem cell does not go on to divide at all. Thus those
% clones will never go extinct. At the same time, (1-theta) of the clones
% behave entirely CP, and follow the usual rules.

basal_bare_average = 1 - 1/gamma*expm1(-gamma*ts);
basal_extinction = generating_function2(r, gamma, ts, 0, 0)';
b = (1-theta) * basal_bare_average ./ (1 - basal_extinction) + ...
    theta * (1 + basal_bare_average);

total_bare_average = basal_bare_average + ...
    (gamma*expm1(-mu*ts) - mu*expm1(-gamma*ts)) / (mu * (mu-gamma));
full_dist = clone_dist_bs_shed(r, gamma, mu, ts, ...
    10, ceil(5 * max(total_bare_average - basal_bare_average)));
[~,n,~] = size(full_dist);
sum_n_p_0n = (0:(n-1)) * squeeze(full_dist(0+1,:,:));
t = (1-theta) * (total_bare_average - sum_n_p_0n) ./ (1 - basal_extinction) + ...
    theta * (total_bare_average + 1);

s = t - b;

end