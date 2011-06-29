function F = posttreat_gf(t1, t2, x0, y0, z0)

r1 = 0.10;
rho1 = 0.65;
lambda1 = 1.87;
gamma1 = rho1/(1-rho1);
mu1 = rho1 / 0.82;

r2 = 0.088;
rho2 = 0.582;
lambda2 = 3.76;
gamma2 = rho2/(1-rho2);
mu2 = rho2 / 1.1;

F = generating_function3_shed(r2, gamma2, mu2, t2*lambda2, ...
        generating_function3_shed(r1, gamma1, mu1, t1*lambda1, x0, y0, z0), ...
        ( (z0-1)*gamma1 * exp(-mu1*lambda1*t1) ...
           + ((y0-z0)*gamma1 + (1-y0)*mu1) * exp(-gamma1*lambda1*t1)) ...
          / (gamma1 - mu1) ...
          + 1, ...
        exp(-mu1*lambda1*t1)*(z0-1) + 1 );

end