function extract_tau(file, taus)

taus = reshape(taus, numel(taus), 1);

% file should contain lambdas, rhos, rs and pxd
load(file);

p = zeros(numel(rs), numel(lambdas), numel(taus));

for i=[1:numel(rs)]
    for j=[1:numel(lambdas)]
        p(i,j,:) = rs(i) * lambdas(j) * interp1(rhos, squeeze(pxd(i,:,j)), taus .* rs(i) .* lambdas(j), [], 0);
    end
end

p_tau = squeeze(trapz(lambdas, trapz(rs, p)));
z = trapz(taus, p_tau);
p_tau = p_tau ./ z;

%newplot(figure);
plot(taus, p_tau);
xlabel('$\tau$/week', 'Interpreter', 'latex');

t1 = trapz(taus, taus .* p_tau);
t2 = trapz(taus, taus.^2 .* p_tau);

fprintf(1, 'mean: %g\tstd: %g\n', t1, sqrt(t2-t1^2));

end