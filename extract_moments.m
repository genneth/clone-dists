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

end
