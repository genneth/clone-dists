function clayton_plot

% The goal is to reproduce the theoretical curves for clone size
% distributions as found in Clayton et al. Nature 2007, Figure 2d. Only
% *survining* clones with more than 1 cell are counted.

nn = 5;
%sz = [2^nn+1 2^nn+1];
sz = [50 50];
P0 = initial_eyp(sz);
lambda = 1.1;
r = 0.08;
rho = 0.22;
gamma = lambda * rho / (1-rho);

[Tl Tr Tg] = generate_transition_matrix(sz);
ts = 10 .^ ([0:40] ./ 20);

pops = arrayfun(@(t)(bin(nn, condPb(Pb(population(Tl, Tr, Tg, lambda, r, gamma, t, P0))))), ts, 'UniformOutput', false);
pops = cell2mat(pops);
gf = newplot;
semilogx(gf, ts, pops);
xlabel(gf, '$t$ / weeks', 'Interpreter', 'latex');
ylabel(gf, 'proportion', 'Interpreter', 'latex');
saveas(gf, 'clayton_plot', 'pdf');

end