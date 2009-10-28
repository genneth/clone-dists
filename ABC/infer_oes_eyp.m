function infer_oes_eyp

t1 = 3.0+1/7; % 8 days
data1 = [2 0 16; 1 1 18; 0 2 20; 3 0 4; 2 1 19; 1 2 10; 0 3 7; 4 0 3; 3 1 8; 2 2 5; 1 3 1; 0 4 4; 4 1 5; 3 2 10; 2 3 3; 1 4 1; 0 5 1; 6 0 1; 5 1 2; 4 2 4; 3 3 3; 2 4 1; 1 5 2; 6 1 1; 5 2 2; 3 4 1; 6 2 3; 5 3 1; 4 4 4; 6 3 1; 5 4 1; 4 5 2; 6 4 1; 9 2 1; 8 3 2; 10 3 1; 9 4 1; 8 5 1; 10 4 1; 11 4 1; 12 4 1; 16 6 1; 18 5 1]';

k = 23;
[Tl Tr Tg] = generate_transition_matrix(k);
P0 = initial_eyp(k);

% x by y by z grid for the parameter space
x = 30;
y = 30;
z = 0;
pdx = zeros(y+1, x+1, z+1);

for i = 1:x+1
    for j = 1:y+1
        for h = 1:z+1
            rho    = (i-1)/x * 0.8 + 0.1;           % 0.1 - 0.9
            r      = (j-1)/y * 0.5;                 % 0.0 - 0.5
            %lambda = (h-1)/z * (1/3 - 1/6) + (1/6); % 3 to 6 week cycle
            lambda = 0.11*7;
            pdx(j,i,h) = PDX(condPbs(Pbs(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t1, P0))), data1);
            %pdx(j,i,h) = pdx(j,i,h) * PDX(condPbs(Pbs(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t3, P0))), data3);
        end
    end
end

rhos = [0:x] / x * 0.8 + 0.1;
rs = [0:y] / y * 0.5;
%lambdas = [0:z] / z * (1/3 - 1/6) + (1/6);
pmax = max(max(max(pdx)));
%isosurface(rhos, rs, lambdas, pdx ./ pmax, 0.5);
contour(rhos, rs, pdx ./ pmax);
xlabel('$\rho$', 'Interpreter', 'latex');
ylabel('$r$', 'Interpreter', 'latex');
%zlabel('$\lambda$', 'Interpreter', 'latex');

end
