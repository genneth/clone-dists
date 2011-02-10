function pbs = clone_dist_bs_shed(r, gamma, mu, ts, M0, N0)

M = M0+1; N = N0+1;

    function pmn = fourier2(fyz)
        pmn = abs(fft(fft(fyz, [], 1), [], 2)) / (M*N);
    end

%curr_samples = zeros(M,N, numel(ts));
curr_ps      = zeros(M,N, numel(ts));
next_samples = zeros(M,N, numel(ts));
next_ps      = zeros(M,N, numel(ts));
for i = 1:M
    parfor j = 1:N
        y = exp(2i*pi/M * (i-1));
        z = exp(2i*pi/N * (j-1));
        next_samples(i,j,:) = generating_function3_shed(r, gamma, mu, ts, y, y, z);
    end
end
next_ps = fourier2(next_samples);

    function e = rel_err
        e = max(max(max(abs(curr_ps(1:(M0+1),1:(N0+1),:) - next_ps(1:(M0+1),1:(N0+1),:)) ./ next_ps(1:(M0+1),1:(N0+1),:))));
    end

    function e = abs_err
        e = max(max(max(abs(curr_ps(1:(M0+1),1:(N0+1),:) - next_ps(1:(M0+1),1:(N0+1),:)))));
    end

while rel_err > 1e-2 && abs_err > 1e-3
    curr_samples = next_samples;
    curr_ps      = next_ps;
    next_s01      = zeros(M,N, numel(ts));
    next_s10      = zeros(M,N, numel(ts));
    next_s11      = zeros(M,N, numel(ts));
    dy = exp(2i*pi/(2*M));
    dz = exp(2i*pi/(2*M));
    for i = 1:M
        parfor j = 1:N
            y = exp(2i*pi/(2*M) * (2*i-2));
            z = exp(2i*pi/(2*N) * (2*j-2));
            next_s01(i,j,:) = generating_function3_shed(r, gamma, mu, ts, y, y, z*dz);
            next_s10(i,j,:) = generating_function3_shed(r, gamma, mu, ts, y*dy, y*dy, z);
            next_s11(i,j,:) = generating_function3_shed(r, gamma, mu, ts, y*dy, y*dy, z*dz);
        end
    end
    % interleave
    next_samples = zeros(2*M,2*N, numel(ts));
    for i = 1:M
        for j = 1:N
            next_samples(2*i-1,2*j-1,:) = curr_samples(i,j,:);
            next_samples(2*i-1,2*j,:)   = next_s01(i,j,:);
            next_samples(2*i,2*j-1,:)   = next_s10(i,j,:);
            next_samples(2*i,2*j,:)     = next_s11(i,j,:);
        end
    end
    M = M*2; N = N*2;
    next_ps = fourier2(next_samples);
%     fprintf(1, 'oversampled %d rel error %g\t abs error %g\n', (M/(M0+1))*(N/(N0+1)), rel_err, abs_err);
end

pbs = next_ps(1:(M0+1), 1:(N0+1), :);

end
