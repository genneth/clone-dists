function pops = exact_pops(lambda, r, gamma, t, N0)

% The key is in choosing the right number of evaluations N. We do this by
% recursively subdividing the unit circle, and watching the convergence.

n = numel(t);
N = N0; % start from the required number
curr_samples = zeros(n, N);
curr_pops = zeros(n, N);
next_samples = zeros(n, N);
next_samples(:,1) = 1;
temp = zeros(n, N-1);
parfor j = [1:(N-1)]
    z = exp(i*2*pi/N * j);
    temp(:,j) = xi(lambda, r, gamma, t, z);
end
for j = [1:(N-1)]
    next_samples(:,j+1) = temp(:,j);
end
next_pops = real(fft(next_samples, [], 2)) ./ N;

    function e = rel_err
        e = max(max(abs(curr_pops(:,1:N0) - next_pops(:,1:N0)) ./ next_pops(:,1:N0)));
    end

    function e = abs_err
        e = max(max(abs(curr_pops(:,1:N0) - next_pops(:,1:N0))));
    end

fprintf(1, 'oversample %d rel error %g\t abs error %g\n', N/N0, rel_err, abs_err);
while rel_err > 1e-3 && abs_err > 1e-7
    curr_samples = next_samples;
    curr_pops = next_pops;
    next2 = zeros(n, N);
    parfor j = [1:N]
        z = exp(i*2*pi/(2*N) * (2*j-1));
        next2(:,j) = xi(lambda, r, gamma, t, z);
    end
    % interleave
    N = N*2;
    next_samples = zeros(n, N);
    for j = [1:N/2]
        next_samples(:, 2*j-1) = curr_samples(:, j);
        next_samples(:, 2*j)   = next2(:,j);
    end
    next_pops = real(fft(next_samples, [], 2)) ./ N;
    fprintf(1, 'oversample %d rel error %g\t abs error %g\n', N/N0, rel_err, abs_err);
end

pops = next_pops(:,1:N0);

end
