function pops = exact_pops(lambda, r, gamma, t, N0)

% The key is in choosing the right number of evaluations N. We do this by
% recursively subdividing the unit circle, and watching the convergence.

N = N0; % start from the required number
curr_samples = zeros(1, N);
curr_pops = zeros(1, N);
next_samples = [1, arrayfun(@(z)(xi(lambda,r,gamma,t,z)), exp(i*2*pi/N.*[1:(N-1)]))];
next_pops = real(fft(next_samples)) ./ N;

    function e = rel_err
        e = max(abs(curr_pops(1:N0) - next_pops(1:N0)) ./ next_pops(1:N0));
    end

    function e = abs_err
        e = max(abs(curr_pops(1:N0) - next_pops(1:N0)));
    end

while rel_err > 1e-3 && abs_err > 1e-7
    curr_samples = next_samples;
    curr_pops = next_pops;
    N = N*2;
    next2 = arrayfun(@(z)(xi(lambda,r,gamma,t,z)), exp(i*2*pi/N.*[1:2:(N-1)]));
    next_samples = reshape([curr_samples ; next2], 1, N); % interleave
    next_pops = real(fft(next_samples)) ./ N;
    fprintf(1, 'rel error %g\t abs error %g\n', rel_err, abs_err);
end

fprintf(1, 'oversampled by %g\n', N/N0);

pops = next_pops(1:N0);

end