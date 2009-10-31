function p = PDX(p0, Pb, data)

p = p0;
for d = data
    m = d(1);
    n = d(2);
    p = p * Pb(m+1)^n;
end

end
