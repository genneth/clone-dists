function P = clone_dist(gamma, q, t, n)

% returns a vector P(i+1) which are the p_i

%n = average_size(gamma, q, t) * 5;
Z = exp(2i*pi/n * (0:(n-1)));
F = arrayfun(@(z)(generating_function(gamma, q, t, z, z)), Z);
P = abs(real(fft(F))) / n;

end