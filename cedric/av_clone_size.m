function [b, s, t] = av_clone_size(r, delta, gamma, mu, ts)

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
% this may thus be computed from generating_function2 and associates.
% 
% p_0n can be given by inverse transforming F(0,0,z) for
% generating_function3, and sum_n n p_0n can be brute-forced.
% 
% The basal average can be computed similarly, with a familiar algorithm.

basal_bare_average = (exp(-gamma*ts) .* (exp((gamma+2*delta)*ts) * (1+gamma) + 2*delta - 1)) / (gamma + 2*delta);
basal_extinction = generating_function2_unbalance(r, delta, gamma, ts, 0, 0)';
b = basal_bare_average ./ (1 - basal_extinction); % sum_n_p_0n

total_bare_average = basal_bare_average + ...
    (gamma*(2*delta-1)*(exp(-mu*ts)*(gamma+2*delta) + exp(2*delta*ts)*(mu-gamma) - exp(-gamma*ts)*(2*delta+mu))) / ((gamma+2*delta)*(gamma-mu)*(2*delta+mu));
extincts = inverse_z_transform(@(z) generating_function3_shed_unbalance(r, delta, gamma, mu, ts, 0, 0, z), ...
    ceil(5 * max(total_bare_average - basal_bare_average)), 1e-3, 1e-4);
[n,~] = size(extincts);
sum_n_p_0n = (0:(n-1)) * extincts;
t = (total_bare_average - sum_n_p_0n) ./ (1 - basal_extinction);

s = t - b;

end