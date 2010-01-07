function moments_plot(sz, lambda, r, gamma)

ts = 10 .^ ([-40:40] ./ 20);
pops = arrayfun(@(t)(exact_pops(lambda, r, gamma, t, sz)), ts, 'UniformOutput', false);
for i = [1:numel(pops)]
    p = pops{i};
    ms(i,1) = sum(p);
    ms(i,2) = sum(([1:numel(p)]'-1) .* p);
end

newplot;
subplot(2, 1, 1); semilogx(ts, ms(:,1));
subplot(2, 1, 2); semilogx(ts, ms(:,2), ...
                           ts, 1+(1 - exp(-gamma .* ts))/(gamma/lambda));

end