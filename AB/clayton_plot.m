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

    function Pbin = bin(Pb)
        Pbin = zeros(nn,1);
        for i = [1:nn]
            Pbin(i) = sum(Pb((2^(i-1)+2):(2^i+1)));
        end
    end

pops = arrayfun(@(t)(bin(condPb(Pb(population(Tl, Tr, Tg, lambda, r, gamma, t, P0))))), ts, 'UniformOutput', false);
pops = cell2mat(pops);
semilogx(ts, pops);
xlabel('$t$ / weeks', 'Interpreter', 'latex');
ylabel('proportion', 'Interpreter', 'latex');

end