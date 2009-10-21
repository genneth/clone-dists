function P0 = initial_edu(k, r)

M = tetra(k+1);
P0 = zeros(M, 1);
P0(t2i([2 0 0])) = r;
P0(t2i([1 1 0])) = 1 - 2*r;
P0(t2i([0 2 0])) = r;

end
