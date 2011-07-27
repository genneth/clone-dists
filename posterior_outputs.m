function posterior_outputs(samples, pf)

constf = @(r,g,l) ones(size(r));
rf = @(r,g,l) r;
rhof = @(r,g,l) g./(1+g);
lf = @(r,g,l) log(l);
tf = @(r,g,l) log(g./(1+g)./(r.*l));
gf = @(r,g,l) log(g.*l);
mf = @(r,g,l) log(g./(1+g)./0.82.*l);
sq = @(f) (@(r,g,l) f(r,g,l).^2);

rs = linear_quadrature_integrate3_multiple({constf, rf, rhof, lf, tf, gf, mf, sq(rf), sq(rhof), sq(lf), sq(tf), sq(gf), sq(mf)}, samples, pf);
Z = rs(1);
r1 = rs(2) / Z;
rho1 = rs(3) / Z;
ll1 = rs(4) / Z;
lt1 = rs(5) / Z;
g1 = rs(6) / Z;
m1 = rs(7) / Z;

% fprintf('r = %.3f, ρ = %.3f, 7/λ = %.2f days, 7*τ = %.2f days\n', r1, rho1, 7/exp(ll1), 7*exp(lt1));

r2 = rs(8) / Z;
rho2 = rs(9) / Z;
ll2 = rs(10) / Z;
lt2 = rs(11) / Z;
g2 = rs(12) / Z;
m2 = rs(13) / Z;

fprintf('r = %.3f ± %.3f\n', r1, sqrt(r2 - r1^2));
fprintf('ρ = %.3f ± %.3f\n', rho1, sqrt(rho2 - rho1^2));
fprintf('λ = %.2f / weeks + %.0f%% - %.0f%%\n', exp(ll1), (exp(sqrt(ll2 - ll1^2)) - 1)*100, (1 - exp(-sqrt(ll2 - ll1^2)))*100);
fprintf('τ = %.2f weeks + %.0f%% - %.0f%%\n', exp(lt1), (exp(sqrt(lt2 - lt1^2)) - 1)*100, (1 - exp(-sqrt(lt2 - lt1^2)))*100);
fprintf('γ = %.3f / weeks + %.0f%% - %.0f%%\n', exp(g1), (exp(sqrt(g2-g1^2)) - 1)*100, (1-exp(-sqrt(g2-g1^2)))*100);
fprintf('μ = %.3f / weeks ± %.0f%% - %.0f%%\n', exp(m1), (exp(sqrt(m2-m1^2)) - 1)*100, (1-exp(-sqrt(m2-m1^2)))*100);

end
