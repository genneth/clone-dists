function samples = sample2(rfun, gammafun, lambdafun, ts2, data2, nsamples, nlambdas, output_file)

% lambdafun is unitful; ts2 is unitful. everything else is unitless.
% rfun, gammafun, lambdafun :: (MonadRandom m) => m Double
% ts2 :: [], data2 :: {} -- times and observations of basal sizes
% samples :: [{r, lambda, gamma, [log p2]}]
%     p2 is likelihoods, not posterior distributions
% nsamples -- number of samples in (gamma,rho) space to take
% nlambdas -- number of lambda points to sample per (gamma,rho) point

assert(numel(ts2) == numel(data2), 'inconsistent sizes between ts2 and data2');

% we use the fact that per invocation of clone_dist_* we can extract many
% separate distributions for different lambdas
samples_r      = zeros(nsamples,nlambdas);
samples_gamma  = zeros(nsamples,nlambdas);
samples_lambda = zeros(nsamples,nlambdas);
samples_p2     = zeros(nsamples,nlambdas,numel(ts2));
samples        = cell(nsamples*nlambdas, 5);

maxN2 = 0;
for i = 1:numel(ts2)
    maxN2 = max(maxN2, numel(data2{i}) - 1);
end

start = tic;
for i = 1:nsamples
    iter_start = tic;
    
    r = feval(rfun);
    lambda  = zeros(nlambdas, 1);
    for j = 1:nlambdas
        lambda(j) = feval(lambdafun);
    end
    gamma = feval(gammafun);

    samples_r(i,:)      = r; 
    samples_gamma(i,:)  = gamma;
    samples_lambda(i,:) = lambda; 

    lts = zeros(numel(ts2) * nlambdas, 1);
    tsj = zeros(numel(ts2) * nlambdas, 1);
    tsk = zeros(numel(ts2) * nlambdas, 1);
    for j = 1:numel(ts2)
        for k = 1:nlambdas
            lts(nlambdas * (j-1) + k) = ts2(j) * lambda(k);
            tsj(nlambdas * (j-1) + k) = j;
            tsk(nlambdas * (j-1) + k) = k;
        end
    end
    [lts, SI] = sort(lts);
    tsj = tsj(SI); tsk = tsk(SI);

    dist2_ = clone_dist_b(r, gamma, lts, maxN2);
    % condition on b>1
    [m_,~] = size(dist2_);
    dist2_ = dist2_ ./ repmat(1 - dist2_(0+1,:) - dist2_(1+1,:), m_, 1);
    dist2_(0+1,:) = 0; dist2_(1+1,:) = 0;

    dist2 = zeros(m_,nlambdas,numel(ts2));
    for l = 1:(numel(ts2)*nlambdas)
        dist2(:,tsk(l),tsj(l)) = dist2_(:,l);
    end

    samples_p2_ = zeros(nlambdas, numel(ts2));
    for j = 1:nlambdas
        for k = 1:numel(ts2)
            [~,cols,v] = find(data2{k});
            for l = 1:numel(v)
                samples_p2_(j,k) = samples_p2_(j,k) + ...
                    v(l) * log(dist2(cols(l),j,k));
            end
        end
    end
    samples_p2(i,:,:) = samples_p2_;
    
    for j = 1:nlambdas
        samples{(i-1)*nlambdas+j,1} = samples_r(i,j);
        samples{(i-1)*nlambdas+j,2} = samples_gamma(i,j);
        samples{(i-1)*nlambdas+j,3} = samples_lambda(i,j);
        samples{(i-1)*nlambdas+j,4} = squeeze(samples_p2(i,j,:));
        samples{(i-1)*nlambdas+j,5} = 1; % placeholder for p3
    end
    partial = samples(1:(i*nlambdas), :);
    save(output_file, 'partial');
    
    fprintf('%d of %d, this iteration: %.1fs, to go: %.2fh\n', i, nsamples, toc(iter_start), (toc(start)/i*(nsamples-i))/3600);
end

save(output_file, 'samples');

end
