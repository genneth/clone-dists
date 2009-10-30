function moments_plot(sz, lambda, r, gamma)

P0 = initial_eyp(sz);
[Tl Tr Tg] = generate_transition_matrix(sz);
ts = 10 .^ ([-40:40] ./ 20);
pops = arrayfun(@(t)(Pb(population(Tl, Tr, Tg, lambda, r, gamma, t, P0))), ts, 'UniformOutput', false);
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