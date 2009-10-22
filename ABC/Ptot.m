function Ptot = Ptot(P)

M = length(P);
k = untetra(M) - 1;
Ptot = zeros(k+1,1);

for i = 1:M
    trip = i2t(i);
    m = trip(1); n = trip(2); l = trip(3);
    Ptot(m+n+l+1) = Ptot(m+n+l+1) + P(i);
end

end
