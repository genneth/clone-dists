function infer_edu

% sanity check with the 2 week EdU data; can be compared with Mathematica

t = 2.0;
data = [1 1 6; 2 0 12; 0 2 2; 2 1 20; 3 0 4; 4 0 1]';

k = 4;
[Tl Tr Tg] = generate_transition_matrix(k);

% n by n grid for the parameter space
n = 30;
pdx = zeros(n+1, n+1);

for i = 1:n+1
    for j = 1:n+1
        rho = (i-1)/n * 0.8 + 0.1; % 0.1 - 0.9
        r   = (j-1)/n * 0.5;       % 0.0 - 0.5
        pdx(j,i) = PDX(Pbs(population(Tl, Tr, Tg, 0.25, r, 0.25 * rho / (1-rho), t, initial_edu(k, r))), data);
    end
end

rhos = [1:n+1] / n * 0.8 + 0.1;
rs = [1:n+1] / n * 0.5;
contour(rhos, rs, pdx);

end
