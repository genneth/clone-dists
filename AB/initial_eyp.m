function P0 = initial_eyp(sz)

P0 = zeros(sz);
P0(p2i(sz,[1 0])) = 1;

end