function posterior_outputs_convergence_test(samples, pf, s)

constf = @(r,g,l) ones(size(r));
rf = @(r,g,l) r;
rhof = @(r,g,l) g./(1+g);
lf = @(r,g,l) log(l);

[n,~] = size(samples);
for i=1:floor(n/s)
    
    rs = linear_quadrature_integrate3_multiple({constf, rf, rhof, lf}, samples(1:(i*s),:), pf);
    Z = rs(1);
    r1(i) = rs(2) / Z;
    rho1(i) = rs(3) / Z;
    ll1(i) = rs(4) / Z;

end

xs = (1:floor(n/s))*s;
plot(xs, r1, xs, rho1, xs, ll1);

end
