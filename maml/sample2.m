function samples = sample2(rfun, gammafun, lambdafun, ts2, data2, nsamples, nlambdas)

% lambdafun is unitful; ts2 is unitful. everything else is unitless.
% rfun, gammafun, lambdafun :: (MonadRandom m) => m (Double, Double)
% ts2 :: [], data2 :: {} -- times and observations of basal sizes
% samples :: [{prob, r, lambda, gamma, [log p2]}]
%     p2 are likelihoods, not posterior distributions
% nsamples -- number of samples in (r,gamma) space to take
% nlambdas -- number of lambda points to sample per (r,gamma) point

assert(numel(ts2) == numel(data2), 'inconsistent sizes between ts2 and data2');

wbh = waitbar(0.0, 'initialising...');

% we use the fact that per invocation of clone_dist_* we can extract many
% separate distributions for different lambdas
samples_prob   = zeros(nsamples,nlambdas);
samples_r      = zeros(nsamples,nlambdas);
samples_gamma  = zeros(nsamples,nlambdas);
samples_lambda = zeros(nsamples,nlambdas);
samples_p2     = zeros(nsamples,nlambdas,numel(ts2));

maxN = 0;
for i = 1:numel(ts2)
    maxN = max(maxN, numel(data2{i}) - 1);
end

waitbar(0.0, wbh, 'running...');
tstart = tic;

for i = 1:nsamples
    
    [rp, r] = feval(rfun);
    lambdap = zeros(nlambdas, 1);
    lambda  = zeros(nlambdas, 1);
    for j = 1:nlambdas
        [lambdap(j), lambda(j)] = feval(lambdafun);
    end
    [gammap, gamma] = feval(gammafun);

    samples_prob(i,:)   = rp*gammap .* lambdap;
    samples_r(i,:)      = r; 
    samples_gamma(i,:)  = gamma;
    samples_lambda(i,:) = lambda; 

%     fprintf(1, 'parameters: r = %f, gamma = %f\n', r, gamma);
%     fprintf(1, '            lambdas =\n');
%     fprintf(1, '                %3f\n', lambda);

    % generate a long list of "times", and remember which time point and
    % lambda it corresponds to. use these "times" to generate the clone
    % size distribution, and then restore to the correct time point and
    % corresponding lambda
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

    dist2_ = clone_dist_b(r, gamma, lts, maxN);
    [m_,~] = size(dist2_);

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
    
    waitbar(i/nsamples, wbh, ...
        sprintf('%.1f%% complete, %.1f min to go', ...
        i/nsamples*100, (toc(tstart)/i*nsamples - toc(tstart))/60));

end

close(wbh);

samples = cell(nsamples*nlambdas, 5);
for i = 1:nsamples
    for j = 1:nlambdas
        samples{(i-1)*nlambdas+j,1} = samples_prob(i,j);
        samples{(i-1)*nlambdas+j,2} = samples_r(i,j);
        samples{(i-1)*nlambdas+j,3} = samples_gamma(i,j);
        samples{(i-1)*nlambdas+j,4} = samples_lambda(i,j);
        samples{(i-1)*nlambdas+j,5} = squeeze(samples_p2(i,j,:));
    end
end

end
