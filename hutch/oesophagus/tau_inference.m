function tau_inference

load intermouse_basal;
epsilon = 0;
ts = [3/7+epsilon*(0:3) 1 10/7+epsilon*(0:3) 3+epsilon*(0:2) 6+epsilon*(0:2) 13+epsilon*(0:3) 26+epsilon*(0:2) 52+epsilon*(0:2)];
t1s = zeros(1, numel(ts));
t2s = zeros(1, numel(ts));
for i=1:numel(ts)
    [t1s(i), t2s(i)] = posterior_output_tau(samples, @(p2s, p3s) p2s(i));
end
dts = sqrt(t2s - t1s.^2);

[t1, t2] = posterior_output_tau(samples, @(p2s, p3s) p2s(i));

fh = figure;
ah = newplot(fh);

errorbar(ah, ts, t1s, dts);


end