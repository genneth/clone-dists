function population_plot(sz, lambda, r, gamma)

P0 = initial_eyp(sz);
[Tl Tr Tg] = generate_transition_matrix(sz);
ts = 10 .^ ([-40:40] ./ 20);
pops = arrayfun(@(t)(Pb(population(Tl, Tr, Tg, lambda, r, gamma, t, P0))), ts, 'UniformOutput', false);
pops = cell2mat(pops);
semilogx(ts, pops);
xlabel('$t$ / weeks', 'Interpreter', 'latex');
ylabel('proportion', 'Interpreter', 'latex');

end