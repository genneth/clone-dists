function P0 = initial_eyp(k)

M = tetra(k+1);
P0 = zeros(M, 1);
P0(t2i([1 0 0])) = 1;

end
