function trip = untriple(z)

w = untetra(z);
t = untriangle(z - tetra(w));
m = z - triangle(t) - tetra(w);
n = t - m;
l = w - t;

trip = [m n l];

end
