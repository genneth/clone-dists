function b = valid_pair(sz, p)

m = p(1); n = p(2);
b = (m >= 0) && (n >= 0) && (m < sz(1)) && (n < sz(2));

end