function Pmnl = Pmnl(P)

M = length(P);
k = untetra(M) - 1;
Pmnl = zeros(k+1,k+1,k+1);

for i = 1:M
    trip = i2t(i);
    m = trip(1); n = trip(2); l = trip(3);
    Pmnl(m+1,n+1,l+1) = P(i);
end

end
