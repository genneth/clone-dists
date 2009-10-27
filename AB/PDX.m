function p = PDX(Pb, data)

p = 1;
for d = data
    m = d(1);
    n = d(2);
    p = p * Pb(m+1)^n;
end

end
