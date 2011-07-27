function result = linear_quadrature_integrate3(integrand, samples, pf)

results = linear_quadrature_integrate3_multiple({integrand}, samples, pf);
result = results(1);

end
