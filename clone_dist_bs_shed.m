function pbs = clone_dist_bs_shed(r, gamma, mu, ts, M0, N0)

pbs = inverse_z_transform2(@(x,y) generating_function3_shed(r, gamma, mu, ts, x, x, y), M0, N0, 1e-3, 1e-4);

end
