function p1 = single_cell(r, gamma, mu, ts)

n = 40;

Z = exp(2i*pi/n * (0:(n-1)));
F0 = zeros(n*n, numel(ts));
parfor i = 0:(n*n-1)
    y = Z(floor(i/n) + 1);
    z = Z(mod(i,n)+1);
    F0(i+1,:) = generating_function(r, gamma, mu, ts, y, y, z);
end
F = zeros(n,n, numel(ts));
for i = 1:n
    for j = 1:n
        F(i,j,:) = F0((i-1)*n + j,:);
    end
end
P = abs(real(fft(fft(F,[],1),[],2))) / n^2;

p1 = squeeze(P(1+1,0+1,:) ./ (1 - P(0+1,0+1,:)));

end
