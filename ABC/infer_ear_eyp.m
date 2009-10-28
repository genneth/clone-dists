function infer_ear_eyp

t1 = 3.0;
data1 = [2 0 11; 1 1 4; 0 2 6; 3 0 1; 2 1 3; 0 3 1; 4 0 1; 4 1 1; 6 1 1]';

t2 = 6.0;
data2 = [2 0 11; 1 1 16; 0 2 21; 3 0 5; 2 1 12; 1 2 7; 0 3 5; 4 0 3; 3 1 8; 2 2 7; 0 4 2; 5 1 2; 4 2 1; 3 4 1]';

t3 = 6.0;
data3 = [2 0 42; 3 0 12; 4 0 7; 6 0 1; 1 1 46; 2 1 38; 3 1 20; 4 1 7; 5 1 2; 0 2 33; 1 2 7; 2 2 12; 3 2 3; 4 2 4; 5 2 2; 0 3 5; 1 3 2; 3 3 1; 4 3 2; 0 4 3; 3 4 1]';

k = 7;
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
            lambda = 0.25;
            pdx(j,i,h) = PDX(condPbs(Pbs(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t1, P0))), data1);
            pdx(j,i,h) = pdx(j,i,h) * PDX(condPbs(Pbs(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t3, P0))), data3);
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
