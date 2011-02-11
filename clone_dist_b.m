function pb = clone_dist_b(r, gamma, ts, N0)

% The key is in choosing the right number of evaluations N. We do this by
% recursively subdividing the unit circle, and watching the convergence.

N = N0+1;
% curr_samples = zeros(N, numel(ts));
curr_ps      = zeros(N, numel(ts));
next_samples = zeros(N, numel(ts));
next_ps      = zeros(N, numel(ts));
parfor j = 1:N
    z = exp(2i*pi/N * (j-1));
    next_samples(j,:) = generating_function2(r, gamma, ts, z, z);
end
next_ps = abs(fft(next_samples, [], 1)) / N;

    function e = rel_err
        e = max(max(abs(curr_ps(1:(N0+1),:) - next_ps(1:(N0+1),:)) ./ next_ps(1:(N0+1),:)));
    end

    function e = abs_err
        e = max(max(abs(curr_ps(1:(N0+1),:) - next_ps(1:(N0+1),:))));
    end

while rel_err > 1e-5 && abs_err > 1e-7
    curr_samples = next_samples;
    curr_ps      = next_ps;
    next_s2      = zeros(N, numel(ts));
    parfor j = 1:N
        z = exp(2i*pi/(2*N) * (2*j-1));
        next_s2(j,:) = generating_function2(r, gamma, ts, z, z);
    end
    % interleave
    N = N*2;
    next_samples = zeros(N, numel(ts));
    for j = 1:N/2
        next_samples(2*j-1,:) = curr_samples(j,:);
        next_samples(2*j,:)   = next_s2(j,:);
    end
    next_ps = abs(fft(next_samples, [], 1)) / N;
%     fprintf(1, 'sampled %d rel error %g\t abs error %g\n', N, rel_err, abs_err);
end

pb = squeeze(next_ps(1:(N0+1),:));

end
