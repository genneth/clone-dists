function pb = clone_dist_b(r, gamma, ts, N0)

pb = inverse_z_transform(@(z) generating_function2(r, gamma, ts, z, z), N0, 1e-4, 1e-5);

end
