% Given a generating function F(z)=sum_n p_n z^n we find the first
% n=[0..N0] probabilities p_n by a suitable contour integral. We require
% F(z) to be defined on the (close) unit disc |z| <= 1. Where possible
% evaluations of F(z) will be carried out in parallel. F(z) can return an
% array for each evaluation. Evaluation is adaptive, and error control is
% via rel_thres, abs_thres.

function pb = inverse_z_transform(F, N0, rel_thres, abs_thres)

% The key is in choosing the right number of evaluations N. We do this by
% recursively subdividing the unit circle, and watching the convergence.

N = N0+1;

% first, find out how many elements there are per each invocation of F
range_shot = feval(F, 1); % this had better equal one...
sz = size(range_shot);

% curr_samples = zeros(N, numel(ts));
curr_ps      = zeros([N sz]);
next_samples = zeros([N sz]);
next_ps      = zeros([N sz]);
next_samples(1,:) = range_shot;
parfor j = 2:N
    z = exp(2i*pi/N * (j-1));
    next_samples(j,:) = feval(F, z);
end
next_ps = abs(fft(next_samples, [], 1)) / N;

    function e = rel_err
        e = max(max(abs(curr_ps(1:(N0+1),:) - next_ps(1:(N0+1),:)) ./ next_ps(1:(N0+1),:)));
    end

    function e = abs_err
        e = max(max(abs(curr_ps(1:(N0+1),:) - next_ps(1:(N0+1),:))));
    end

while rel_err > rel_thres && abs_err > abs_thres
    curr_samples = next_samples;
    curr_ps      = next_ps;
    next_s2      = zeros([N sz]);
    parfor j = 1:N
        z = exp(2i*pi/(2*N) * (2*j-1));
        next_s2(j,:) = feval(F, z);
    end
    % interleave
    N = N*2;
    next_samples = zeros([N sz]);
    for j = 1:N/2
        next_samples(2*j-1,:) = curr_samples(j,:);
        next_samples(2*j,:)   = next_s2(j,:);
    end
    next_ps = abs(fft(next_samples, [], 1)) / N;
%     fprintf(1, 'sampled %d rel error %g\t abs error %g\n', N, rel_err, abs_err);
end

pb = squeeze(next_ps(1:(N0+1),:));

end
