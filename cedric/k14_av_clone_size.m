function [b, s, t] = k14_av_clone_size(r, delta, gamma, mu, theta, ts)

% For k14, we make the following assumption: theta is the proportion of
% induced clones which are stem like, and they immediatedly (possible with
% some time delay) divide to give another stem cell and a CP cell: S ->
% S+A, and then continue to do so at the correct rate to replace the CP
% loss. Thus those clones will never go extinct. At the same time,
% (1-theta) of the clones behave entirely CP, and follow the usual rules.

cp_basal_bare_average = (exp(-gamma*ts) .* (exp((gamma+2*delta)*ts) * (1+gamma) + 2*delta - 1)) / (gamma + 2*delta);
cp_basal_extinction = generating_function2_unbalance(r, delta, gamma, ts, 0, 0)';
b = (1-theta) * cp_basal_bare_average ./ (1 - cp_basal_extinction) + ...
    theta * (2-((1-exp(-gamma*ts))*(2*delta-1))/gamma);

cp_total_bare_average = cp_basal_bare_average + ...
    (gamma*(2*delta-1)*(exp(-mu*ts)*(gamma+2*delta) + exp(2*delta*ts)*(mu-gamma) - exp(-gamma*ts)*(2*delta+mu))) / ((gamma+2*delta)*(gamma-mu)*(2*delta+mu));
cp_extincts = inverse_z_transform(@(z) generating_function3_shed_unbalance(r, delta, gamma, mu, ts, 0, 0, z), ...
    ceil(5 * max(cp_total_bare_average - cp_basal_bare_average)), 1e-3, 1e-4);
[n,~] = size(cp_extincts);
sum_n_p_0n = (0:(n-1)) * cp_extincts;
t = (1-theta) * (cp_total_bare_average - sum_n_p_0n) ./ (1 - cp_basal_extinction) + ...
    theta * ((2+(1/mu+1/gamma)*(1-2*delta)) + ((2*delta-1)*(gamma^2*exp(-mu*ts) - mu^2*exp(-gamma*ts)))/(gamma*mu*(gamma-mu)));

s = t - b;

end