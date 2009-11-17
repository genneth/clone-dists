function pops = exact_pops(lambda, r, gamma, t, N0)

if N0 < 32
    N = 32;
else
    N = N0;
end

samples = arrayfun(@(z)(xi(lambda,r,gamma,t,z)), exp(i*2*pi/N.*[1:(N-1)]));
pops = real(fft([1,samples])) ./ N;

pops = pops(1:N0);

end