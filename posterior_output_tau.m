function [t, dt] = posterior_output_tau(samples, pf)

constf = @(r,g,l) ones(size(r));
tf = @(r,g,l) log(g./(1+g)./(r.*l));
sq = @(f) (@(r,g,l) f(r,g,l).^2);

rs = linear_quadrature_integrate3_multiple({constf, tf, sq(tf)}, samples, pf);
Z = rs(1);
t = rs(2) / Z;
t2 = rs(3) / Z;
dt = sqrt(t2-t*t);

end
