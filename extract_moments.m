function [rho1,rho2,r1,r2,lambda1,lambda2,gamma1,gamma2] = extract_moments(file)

% file should contain lambdas, rhos, rs and pxd
load(file);

rho1 = trapz(rhos, rhos .* trapz(lambdas, trapz(rs, pxd), 3));
rho2 = trapz(rhos, rhos.^2 .* trapz(lambdas, trapz(rs, pxd), 3));

r1 = trapz(rs, rs' .* trapz(lambdas, trapz(rhos, pxd, 2), 3));
r2 = trapz(rs, rs'.^2 .* trapz(lambdas, trapz(rhos, pxd, 2), 3));

lambda1 = trapz(lambdas, lambdas' .* squeeze(trapz(rhos, trapz(rs, pxd))));
lambda2 = trapz(lambdas, lambdas'.^2 .* squeeze(trapz(rhos, trapz(rs, pxd))));

gamma1 = trapz(lambdas, trapz(rhos, ((rhos ./ (1-rhos))' * lambdas) .* squeeze(trapz(rs, pxd))));
gamma2 = trapz(lambdas, trapz(rhos, ((rhos ./ (1-rhos))' * lambdas).^2 .* squeeze(trapz(rs, pxd))));

fprintf(1, 'parameters (+- 95%% confidence intervals):\n');
fprintf(1, '  1/lambda = %.1f +- %.1f days\n', 7/lambda1, 7/lambda1* 2*sqrt(lambda2-lambda1^2)/lambda1);
fprintf(1, '  1/gamma = %.1f +- %.1f days\n', 7/gamma1, 7/gamma1* 2*sqrt(gamma2-gamma1^2)/gamma1);
fprintf(1, '  rho = %.2f +- %.2f\n', rho1, 2*sqrt(rho2-rho1^2));
fprintf(1, '  r = %.2f +- %.2f\n', r1, 2*sqrt(r2-r1^2));

end
