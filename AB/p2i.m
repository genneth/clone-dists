function i = p2i(sz, p)

p1 = p + [1 1];
i = sub2ind(sz, p1(1), p1(2));

end