function p = i2p(sz, i)

[m1,n1] = ind2sub(sz, i);
p = [m1-1 n1-1];

end