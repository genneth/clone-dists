function posterior_outputs(samples, pf)

constf = @(r,g,l) 1;
rf = @(r,g,l) r;
rhof = @(r,g,l) g/(1+g);
lf = @(r,g,l) log(l);
tf = @(r,g,l) log(g/(1+g)/(r*l));
gf = @(r,g,l) log(g*l);
mf = @(r,g,l) log(g/(1+g)/0.82*l);
sq = @(f) (@(r,g,l) f(r,g,l)^2);

Z = linear_quadrature_integrate3(constf, samples, pf);
r1 = linear_quadrature_integrate3(rf, samples, pf) / Z;
rho1 = linear_quadrature_integrate3(rhof, samples, pf) / Z;
ll1 = linear_quadrature_integrate3(lf, samples, pf) / Z;
lt1 = linear_quadrature_integrate3(tf, samples, pf) / Z;
g1 = linear_quadrature_integrate3(gf, samples, pf) / Z;
m1 = linear_quadrature_integrate3(mf, samples, pf) / Z;

% fprintf('r = %.3f, ρ = %.3f, 7/λ = %.2f days, 7*τ = %.2f days\n', r1, rho1, 7/exp(ll1), 7*exp(lt1));

r2 = linear_quadrature_integrate3(sq(rf), samples, pf) / Z;
rho2 = linear_quadrature_integrate3(sq(rhof), samples, pf) / Z;
ll2 = linear_quadrature_integrate3(sq(lf), samples, pf) / Z;
lt2 = linear_quadrature_integrate3(sq(tf), samples, pf) / Z;
g2 = linear_quadrature_integrate3(sq(gf), samples, pf) / Z;
m2 = linear_quadrature_integrate3(sq(mf), samples, pf) / Z;

fprintf('r = %.3f ± %.3f\n', r1, sqrt(r2 - r1^2));
fprintf('ρ = %.3f ± %.3f\n', rho1, sqrt(rho2 - rho1^2));
fprintf('λ = %.2f / weeks + %.0f%% - %.0f%%\n', exp(ll1), (exp(sqrt(ll2 - ll1^2)) - 1)*100, (1 - exp(-sqrt(ll2 - ll1^2)))*100);
fprintf('τ = %.2f weeks + %.0f%% - %.0f%%\n', exp(lt1), (exp(sqrt(lt2 - lt1^2)) - 1)*100, (1 - exp(-sqrt(lt2 - lt1^2)))*100);
fprintf('γ = %.3f / weeks + %.0f%% - %.0f%%\n', exp(g1), (exp(sqrt(g2-g1^2)) - 1)*100, (1-exp(-sqrt(g2-g1^2)))*100);
fprintf('μ = %.3f / weeks ± %.0f%% - %.0f%%\n', exp(m1), (exp(sqrt(m2-m1^2)) - 1)*100, (1-exp(-sqrt(m2-m1^2)))*100);

end