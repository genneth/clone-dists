function P = clone_dist(r1, r2, ts, n)

Z = exp(2i*pi/n * (0:(n-1)));
F = zeros(n, numel(ts));
parfor i = 0:(n-1)
    F(i+1,:) = generating_function(r1, r2, ts, Z(i+1), Z(i+1));
end
P = abs(real(fft(F))) / n;

end
