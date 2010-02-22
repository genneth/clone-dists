function apply_prior(file_stem)

load(strcat(file_stem, '.mat'));

pxd = repmat(reshape(1 ./ lambdas, 1, 1, numel(lambdas)), numel(rs), numel(rhos)) .* pxd;
Z = trapz(lambdas, trapz(rhos, trapz(rs, pxd)));
pxd = pxd ./ Z;

save(strcat(file_stem, '_prior.mat'), 'rhos', 'rs', 'lambdas', 'pxd');

[rho1,rho2,r1,r2,lambda1,lambda2] = extract_moments(strcat(file_stem, '_prior.mat'));
fprintf('rho:\t%g +- %g\n', rho1, sqrt(rho2 - rho1^2));
fprintf('r:\t%g +- %g\n', r1, sqrt(r2 - r1^2));
fprintf('lambda:\t%g +- %g\n', lambda1, sqrt(lambda2 - lambda1^2));

end