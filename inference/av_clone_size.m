function [b, s, t] = av_clone_size(r, gamma, mu, ts)

% The funny thing about averages is that there are so many. We want the
% average basal, suprabasal and total clone sizes, but with a small caveat:
% we only want to count clones which are "rooted", i.e. have a basal count
% of at least one.
% 
% For the totals, we have: 
%   "<m+n>" = sum_{m>0,n} (m+n) p_mn / (1 - sum_n p_0n)
%           = (<m+n> - sum_n n p_0n) / (1 - sum_n p_0n)
% The bare average <m+n> is simple to evaluate, and has a closed form. The
% two summations are more problematic. The latter, sum_n p_0n may be
% seen as just the extinction probability for the basal-layer-only process;
% this may thus be computed from generating_function2.
%
% Currently, we use the brute-force approach to sum_n n p_0n, by evaluating
% the full p_mn and computing the sum directly.
% 
% The basal average can be computed similarly, with a familiar algorithm.

basal_bare_average = 1 - 1/gamma*expm1(-gamma*ts);
basal_extinction = generating_function2(r, gamma, ts, 0, 0)';
b = basal_bare_average ./ (1 - basal_extinction); % sum_n_p_0n

total_bare_average = basal_bare_average + ...
    (gamma*expm1(-mu*ts) - mu*expm1(-gamma*ts)) / (mu * (mu-gamma));
full_dist = clone_dist_bs_shed(r, gamma, mu, ts, ...
    10, ceil(5 * max(total_bare_average - basal_bare_average)));
[~,n,~] = size(full_dist);
sum_n_p_0n = (0:(n-1)) * squeeze(full_dist(0+1,:,:));
t = (total_bare_average - sum_n_p_0n) ./ (1 - basal_extinction);

s = t - b;

end