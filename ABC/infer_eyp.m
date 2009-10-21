function infer_eyp

%t = 3.0;
%data = [2 0 11; 1 1 4; 0 2 6; 3 0 1; 2 1 3; 0 3 1; 4 0 1; 4 1 1; 6 1 1]';

t = 6.0;
data = [2 0 11; 1 1 16; 0 2 21; 3 0 5; 2 1 12; 1 2 7; 0 3 5; 4 0 3; 3 1 8; 2 2 7; 0 4 2; 5 1 2; 4 2 1; 3 4 1]';

k = 7;
[Tl Tr Tg] = generate_transition_matrix(k);
P0 = initial_eyp(k);

% n by n grid for the parameter space
n = 30;
pdx = zeros(n+1, n+1);

for i = 1:n+1
    for j = 1:n+1
        rho = (i-1)/n * 0.8 + 0.1; % 0.1 - 0.9
        r   = (j-1)/n * 0.5;       % 0.0 - 0.5
        pdx(j,i) = PDX(condPbs(Pbs(population(Tl, Tr, Tg, 0.25, r, 0.25 * rho / (1-rho), t, P0))), data);
    end
end

rhos = [1:n+1] / n * 0.8 + 0.1;
rs = [1:n+1] / n * 0.5;
contour(rhos, rs, pdx);

end
