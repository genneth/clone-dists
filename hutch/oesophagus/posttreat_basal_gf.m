function F = posttreat_basal_gf(t1, t2, x0, y0)

r1 = 0.10;
rho1 = 0.65;
lambda1 = 1.87;
gamma1 = rho1/(1-rho1);

r2 = 0.088;
rho2 = 0.582;
lambda2 = 3.76;
gamma2 = rho2/(1-rho2);

F = generating_function2(r2, gamma2, t2*lambda2, ...
        generating_function2(r1, gamma1, t1*lambda1, x0, y0), ...
        exp(-gamma1*lambda1*t1)*(y0-1) + 1);

end