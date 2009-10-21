function Pbs = Pbs(P)

M = length(P);
k = untetra(M) - 1;
Pbs = zeros(k+1,k+1);

for i = 1:M
    trip = i2t(i);
    m = trip(1); n = trip(2); l = trip(3);
    Pbs(m+n+1, l+1) = Pbs(m+n+1, l+1) + P(i);
end

end
