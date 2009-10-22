function clayton_plot

% The goal is to reproduce the theoretical curves for clone size
% distributions as found in Clayton et al. Nature 2007, Figure 2b. Only
% *survining* clones with more than 1 cell are counted. In the current 
% model, that means clones with >0 basal cells and >1 cells overall

nn = 5;
k = 2^nn; % need a power of 2
P0 = initial_eyp(k); % one A cell to begin with
lambda = 1.1;
r = 0.08;
rho = 0.22;
gamma = lambda * rho / (1-rho);

[Tl Tr Tg] = generate_transition_matrix(k);
ts = 10 .^ ([0:20] ./ 20);

    function Pmnl2 = normalise(Pmnl)
        extinct = sum(Pmnl(0+1,0+1,:));
        single_basal = Pmnl(1+1,0+1,0+1) + Pmnl(0+1,1+1,0+1);
        Pmnl(0+1,0+1,:) = 0;
        Pmnl(1+1,0+1,0+1) = 0;
        Pmnl(0+1,1+1,0+1) = 0;
        Pmnl2 = Pmnl ./ (1 - extinct - single_basal);
    end

    function Ptot = collapse(Pmnl)
        ps = size(Pmnl);
        Ptot = zeros(sum(ps)-2, 1);
        for m = 0:(ps(1)-1)
            for n = 0:(ps(2)-1)
                for l = 0:(ps(3)-1)
                    Ptot(m+n+l+1) = Ptot(m+n+l+1) + Pmnl(m+1,n+1,l+1);
                end
            end
        end
    end

pops = arrayfun(@(t)(collapse(normalise(Pmnl(population(Tl, Tr, Tg, lambda, r, gamma, t, P0))))), ts, 'UniformOutput', false);
pops = cell2mat(pops)';
spop = zeros(length(ts), nn);
for i = [1:nn]
    spop(:,i) = sum(pops(:,(2^(i-1)+2):(2^i+1)),2);
end
semilogx(ts, spop);

end
