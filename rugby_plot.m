function rugby_plot(gf, file)

load(file);

[rho1,rho2, r1,r2, lambda1,lambda2] = extract_moments(file);

contourslice(gf, rhos,rs,lambdas, max(0,pxd), rho1,r1,lambda1);
xlabel('$\rho$', 'Interpreter', 'latex');
ylabel('$r$', 'Interpreter', 'latex');
zlabel('$\lambda$', 'Interpreter', 'latex');

end
