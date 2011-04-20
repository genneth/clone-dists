function posterior_outputs(samples, pf)

constf = @(r,g,l) 1;
rf = @(r,g,l) r;
rhof = @(r,g,l) g/(1+g);
lf = @(r,g,l) l;
tf = @(r,g,l) g/(1+g)/(r*l);
sq = @(f) (@(r,g,l) f(r,g,l)^2);

Z = linear_quadrature_integrate3(constf, samples, pf);
r1 = linear_quadrature_integrate3(rf, samples, pf) / Z;
rho1 = linear_quadrature_integrate3(rhof, samples, pf) / Z;
l1 = linear_quadrature_integrate3(lf, samples, pf) / Z;
t1 = linear_quadrature_integrate3(tf, samples, pf) / Z;

fprintf('r = %.3f, ρ = %.3f, 7/λ = %.2f days, 7*τ = %.2f days\n', r1, rho1, 7/l1, 7*t1);

r2 = linear_quadrature_integrate3(sq(rf), samples, pf) / Z;
rho2 = linear_quadrature_integrate3(sq(rhof), samples, pf) / Z;
l2 = linear_quadrature_integrate3(sq(lf), samples, pf) / Z;
t2 = linear_quadrature_integrate3(sq(tf), samples, pf) / Z;

fprintf('r = %.3f ± %.3f\n', r1, sqrt(r2 - r1^2));
fprintf('ρ = %.3f ± %.3f\n', rho1, sqrt(rho2 - rho1^2));
fprintf('λ = %.2f ± %.2f / week\n', l1, sqrt(l2 - l1^2));
fprintf('τ = %.2f ± %.2f weeks\n', t1, sqrt(t2 - t1^2));

end