function k = triple(trip)

m = trip(1);
n = trip(2);
l = trip(3);
k = tetra(m+n+l) + triangle(m+n) + m;

end
