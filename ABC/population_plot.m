function population_plot(k, lambda, r, gamma)

P0 = initial_eyp(k);
[Tl Tr Tg] = generate_transition_matrix(k);
pops = arrayfun(@(x)(population(Tl, Tr, Tg, lambda, r, gamma, 10^(x/10), P0)),[-20:20], 'UniformOutput', false);
pops = cell2mat(pops)';
plot(pops);

end
