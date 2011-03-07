% Given a generating function F(x,y)=sum_m,n p_mn x^m y^n we find the
% probabilities p_mn (0<=m<=M0, 0<=n<=N0) by a suitable contour integral.
% We require F(x,y) to be defined on the (close) volume |x| <= 1 and |y| <=
% 1. Where possible evaluations of F will be carried out in parallel. F can
% return an array for each evaluation. Evaluation is adaptive, and error
% control is via rel_thres, abs_thres.

function pbs = inverse_z_transform2(F, M0, N0, rel_thres, abs_thres)

M = M0+1; N = N0+1;

    function pmn = fourier2(fyz)
        pmn = abs(fft(fft(fyz, [], 1), [], 2)) / (M*N);
    end

% first, find out how many elements there are per each invocation of F
range_shot = feval(F, 1, 1); % this had better equal one...
sz = size(range_shot);

%curr_samples = zeros(M,N, numel(ts));
curr_ps      = zeros([M N sz]);
next_samples = zeros([M N sz]);
next_ps      = zeros([M N sz]);
for i = 1:M
    parfor j = 1:N
        y = exp(2i*pi/M * (i-1));
        z = exp(2i*pi/N * (j-1));
        next_samples(i,j,:) = feval(F, y, z);
    end
end
next_ps = fourier2(next_samples);

    function e = rel_err
        e = max(max(max(abs(curr_ps(1:(M0+1),1:(N0+1),:) - next_ps(1:(M0+1),1:(N0+1),:)) ./ next_ps(1:(M0+1),1:(N0+1),:))));
    end

    function e = abs_err
        e = max(max(max(abs(curr_ps(1:(M0+1),1:(N0+1),:) - next_ps(1:(M0+1),1:(N0+1),:)))));
    end

while rel_err > rel_thres && abs_err > abs_thres
    curr_samples = next_samples;
    curr_ps      = next_ps;
    next_s01      = zeros([M N sz]);
    next_s10      = zeros([M N sz]);
    next_s11      = zeros([M N sz]);
    dy = exp(2i*pi/(2*M));
    dz = exp(2i*pi/(2*M));
    for i = 1:M
        parfor j = 1:N
            y = exp(2i*pi/(2*M) * (2*i-2));
            z = exp(2i*pi/(2*N) * (2*j-2));
            next_s01(i,j,:) = feval(F, y, z*dz);
            next_s10(i,j,:) = feval(F, y*dy, z);
            next_s11(i,j,:) = feval(F, y*dy, z*dz);
        end
    end
    % interleave
    next_samples = zeros([2*M 2*N sz]);
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

pbs = abs(next_ps(1:(M0+1), 1:(N0+1), :));

end
