function samples = sample3_no_shedding(rfun, gammafun, lambdafun, ts2, data2, ts3, data3, nsamples, nlambdas)

% lambdafun is unitful; ts2 and ts3 are unitful. everything else is unitless.
% rfun, gammafun, lambdafun :: (MonadRandom m) => m (Double, Double)
% ts2 :: [], data2 :: {} -- times and observations of basal sizes
% ts3 :: [], data3 :: {} -- times and observations of suprabasal clone sizes
% samples :: [{prob, r, lambda, gamma, [log p2], [log p3]}]
%     p2 and p3 are likelihoods, not posterior distributions
% nsamples -- number of samples in (gamma,rho) space to take
% nlambdas -- number of lambda points to sample per (gamma,rho) point

assert(numel(ts2) == numel(data2), 'inconsistent sizes between ts2 and data2');
assert(numel(ts3) == numel(data3), 'inconsistent sizes between ts3 and data3');
assert(numel(ts3) > 0, 'no suprabasal data. should run sample2');

% wbh = waitbar(0.0, 'initialising...');

% we use the fact that per invocation of clone_dist_* we can extract many
% separate distributions for different lambdas
samples_prob   = zeros(nsamples,nlambdas);
samples_r      = zeros(nsamples,nlambdas);
samples_gamma  = zeros(nsamples,nlambdas);
samples_lambda = zeros(nsamples,nlambdas);
if numel(ts2) > 0
    samples_p2     = zeros(nsamples,nlambdas,numel(ts2));
else
    samples_p2     = zeros(nsamples,nlambdas,1);
end
samples_p3     = zeros(nsamples,nlambdas,numel(ts3));

maxN2 = 0;
for i = 1:numel(ts2)
    maxN2 = max(maxN2, numel(data2{i}) - 1);
end

maxK = 0;
for i = 1:numel(ts3)
    [rows,cols,~] = find(data3{i});
    maxK = max(maxK, max(rows + cols) - 2);
end

[Tl Tr Tg P0] = clone_dist_bs_expv_setup(maxK);

% waitbar(0.0, wbh, 'running...');
% tstart = tic;

parfor i = 1:nsamples
    
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
    lts = zeros(numel(ts3) * nlambdas, 1);
    tsj = zeros(numel(ts3) * nlambdas, 1);
    tsk = zeros(numel(ts3) * nlambdas, 1);
    for j = 1:numel(ts3)
        for k = 1:nlambdas
            lts(nlambdas * (j-1) + k) = ts3(j) * lambda(k);
            tsj(nlambdas * (j-1) + k) = j;
            tsk(nlambdas * (j-1) + k) = k;
        end
    end
    [lts, SI] = sort(lts);
    tsj = tsj(SI); tsk = tsk(SI);
    
    dist3_ = clone_dist_bs_expv_cond(Tl, Tr, Tg, P0, r, gamma, lts);
    [m_,n_,~] = size(dist3_);
    
    dist3 = zeros(m_,n_,nlambdas,numel(ts3));
    for l = 1:(numel(ts3)*nlambdas)
        dist3(:,:,tsk(l),tsj(l)) = dist3_(:,:,l);
    end

    samples_p3_ = zeros(nlambdas,numel(ts3));
    for j = 1:nlambdas
        for k = 1:numel(ts3)
            [row,col,v] = find(data3{k});
            for l = numel(v)
                samples_p3_(j,k) = ...
                    samples_p3_(j,k) ...
                  + log(dist3(row(l),col(l),j,k)) * data3{k}(row(l),col(l));
            end
        end
    end
    samples_p3(i,:,:) = samples_p3_;
    
    % and again, with the basal data, if we have it
    if numel(ts2) > 0
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
                [~,col,v] = find(data2{k});
                for l = numel(v)
                    samples_p2_(j,k) = ...
                        samples_p2_(j,k) ...
                      + log(dist2(col(l),j,k)) * data2{k}(col(l));
                end
            end
        end
    else
        samples_p2_ = ones(mult, 1);
    end
    samples_p2(i,:,:) = samples_p2_;
    
    
%     waitbar(i/nsamples, wbh, ...
%         sprintf('%.1f%% complete, %.1f min to go', ...
%         i/nsamples*100, (toc(tstart)/i*nsamples - toc(tstart))/60));

end

% close(wbh);

samples = cell(nsamples*nlambdas, 6);
for i = 1:nsamples
    for j = 1:nlambdas
        samples{(i-1)*nlambdas+j,1} = samples_prob(i,j);
        samples{(i-1)*nlambdas+j,2} = samples_r(i,j);
        samples{(i-1)*nlambdas+j,3} = samples_gamma(i,j);
        samples{(i-1)*nlambdas+j,4} = samples_lambda(i,j);
        samples{(i-1)*nlambdas+j,5} = squeeze(samples_p2(i,j,:));
        samples{(i-1)*nlambdas+j,6} = squeeze(samples_p3(i,j,:));
    end
end

end
