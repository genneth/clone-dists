function [t1, t2] = posterior_output_tau(samples, pf)

constf = @(r,g,l) 1;
tf = @(r,g,l) log(g/(1+g)/(r*l));
sq = @(f) (@(r,g,l) f(r,g,l)^2);

Z = linear_quadrature_integrate3(constf, samples, pf);
t1 = linear_quadrature_integrate3(tf, samples, pf) / Z;
t2 = linear_quadrature_integrate3(sq(tf), samples, pf) / Z;

end