function infer_ear_eyp

% Basal layer data on ear epidermis. Should very sensitive to the scaling
% limit, but be significantly worse for distinguishing transient behaviour.
% Works almost orthogonally to the short term basal+suprabasal experiment.

t1 = 3;
data1 = [2 65; 3 19; 4 8; 5 1; 6 1]';

t2 = 6;
data2 = [2 112; 3 53; 4 27; 5 17; 6 6; 7 1; 8 3]';

t3 = 26;
data3 = [2 57; 3 38; 4 20; 5 24; 6 15; 7 8; 8 2; 9 8; 10 6; 11 2; 12 2; 13 5; 16 1; 20 1; 22 1]';

t4 = 52;
data4 = [2 33; 3 26; 4 21; 5 15; 6 10; 7 9; 8 6; 9 8; 10 6; 11 1; 12 5; 13 6; 14 2; 15 3; 16 7; 17 2; 19 1; 20 1; 21 1; 25 2; 35 1; 44 1; 50 1]';

sz = [100 100];
[Tl Tr Tg] = generate_transition_matrix(sz);
P0 = initial_eyp(sz);

% x by y by z grid for the parameter space
x = 30;
y = 30;
z = 0;
pdx = zeros(y+1, x+1, z+1);

for i = 1:x+1
    for j = 1:y+1
        for h = 1:z+1
            rho    = (i-1)/x * 0.8 + 0.1;           % 0.1 - 0.9
            r      = (j-1)/y * 0.5 + 0.0;          % 0.0 - 0.5
            %lambda = (h-1)/z * (1/3 - 1/6) + (1/6); % 3 to 6 week cycle
            lambda = 0.25;
            pdx(j,i,h) = 1e200;
            pdx(j,i,h) = pdx(j,i,h) * PDX(condPb(Pb(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t1, P0))), data1);
            pdx(j,i,h) = pdx(j,i,h) * PDX(condPb(Pb(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t2, P0))), data2);
            pdx(j,i,h) = pdx(j,i,h) * PDX(condPb(Pb(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t3, P0))), data3);
            pdx(j,i,h) = 1e50 * pdx(j,i,h) * PDX(condPb(Pb(population(Tl, Tr, Tg, lambda, r, lambda * rho / (1-rho), t4, P0))), data4);
        end
    end
end

rhos = [0:x] / x * 0.8 + 0.1;
rs = [0:y] / y * 0.5 + 0.0;
%lambdas = [0:z] / z * (1/3 - 1/6) + (1/6);
pmax = max(max(max(pdx)));
%isosurface(rhos, rs, lambdas, pdx ./ pmax, 0.5);
%contourslice(rhos, rs, lambdas, pdx ./ pmax, [0.6], [0.25], [0.25]);
contour(rhos, rs, pdx ./ pmax);
xlabel('$\rho$', 'Interpreter', 'latex');
ylabel('$r$', 'Interpreter', 'latex');
%zlabel('$\lambda$', 'Interpreter', 'latex');

end
