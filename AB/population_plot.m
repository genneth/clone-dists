function population_plot(sz, lambda, r, gamma)

ts = 10 .^ ([-40:40] ./ 20);
pops = arrayfun(@(t)(exact_pops(lambda, r, gamma, t, sz)), ts, 'UniformOutput', false);
pops = cell2mat(pops);
semilogx(ts, pops);
xlabel('$t$ / weeks', 'Interpreter', 'latex');
ylabel('proportion', 'Interpreter', 'latex');

end