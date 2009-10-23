function population_plot(k, lambda, r, gamma)

P0 = initial_eyp(k); % single A cell
[Tl Tr Tg] = generate_transition_matrix(k);
ts = 10 .^ ([-40:40] ./ 20);
pops = arrayfun(@(t)(population(Tl, Tr, Tg, lambda, r, gamma, t, P0)), ts, 'UniformOutput', false);
pops = cell2mat(pops)';
semilogx(ts, pops);
xlabel('$t$ / weeks', 'Interpreter', 'latex');
ylabel('proportion', 'Interpreter', 'latex');

end
