function rugby_plot(gf, file)

load(file);

isosurface(rhos,rs,lambdas,pxd, find_cutoff(pxd, 0.683));
isosurface(rhos,rs,lambdas,pxd, find_cutoff(pxd, 0.95));
alpha(0.4)
% hold on
% pxd = squeeze(trapz(lambdas, pxd, 3));
% contour(gca, rhos, rs, pxd, [find_cutoff(pxd, 0.683) find_cutoff(pxd, 0.95)]);
% hold off
% [rho1,rho2, r1,r2, lambda1,lambda2] = extract_moments(file);
% 
% contourslice(gf, rhos,rs,lambdas, max(0,pxd), rho1,r1,lambda1);
xlabel('$\rho$', 'Interpreter', 'latex');
ylabel('$r$', 'Interpreter', 'latex');
zlabel('$\lambda$', 'Interpreter', 'latex');

end
