function clayton_plot

% The goal is to reproduce the theoretical curves for clone size
% distributions as found in Clayton et al. Nature 2007, Figure 2d. Only
% *survining* clones with more than 1 cell are counted.

nn = 5;
N = 2^nn + 1;
lambda = 1.1;
r = 0.08;
rho = 0.22;
gamma = lambda * rho / (1-rho);

ts = 10 .^ ([0:40] ./ 20);

pops = arrayfun(@(t)(bin(nn, condPb(exact_pops(lambda, r, gamma, t, N)))), ts, 'UniformOutput', false);
pops = cell2mat(pops);
gf = newplot;
semilogx(gf, ts, pops);
xlabel(gf, '$t$ / weeks', 'Interpreter', 'latex');
ylabel(gf, 'proportion', 'Interpreter', 'latex');
saveas(gf, 'clayton_plot', 'pdf');

end