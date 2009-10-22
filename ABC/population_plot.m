function population_plot(k, lambda, r, gamma)

P0 = initial_eyp(k);
[Tl Tr Tg] = generate_transition_matrix(k);
ts = 10 .^ ([0:20] ./ 20);
pops = arrayfun(@(t)(condPtot(Ptot(population(Tl, Tr, Tg, lambda, r, gamma, t, P0)))),ts, 'UniformOutput', false);
pops = cell2mat(pops)';
size(pops)
spop = zeros(length(ts), 4);
spop(:,1) = pops(:,2+1);
spop(:,2) = pops(:,3+1)+pops(:,4+1);
spop(:,3) = pops(:,5+1)+pops(:,6+1)+pops(:,7+1)+pops(:,8+1);
spop(:,4) = pops(:,9+1)+pops(:,10+1)+pops(:,11+1)+pops(:,12+1)+pops(:,13+1)+pops(:,14+1)+pops(:,15+1)+pops(:,16+1);
semilogx(ts, spop);

end
