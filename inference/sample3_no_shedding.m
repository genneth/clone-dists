function samples = sample3_no_shedding(rfun, gammafun, lambdafun, ts2, data2, ts3, data3, n, ppm)

if nargin == 8
    ppm_ = 0;
elseif nargin == 9
    ppm_ = 1;
end

% lambdafun is unitful; ts2 and ts3 are unitful. everything else is unitless.
% rfun, gammafun, lambdafun :: (MonadRandom m) => m (Double, Double)
% ts2 :: [], data2 :: {} -- times and observations of basal sizes
% ts3 :: [], data3 :: {} -- times and observations of suprabasal clone sizes
% samples :: [{prob, r, lambda, gamma, [log p2], [log p3]}]
%     p2 and p3 are likelihoods, not posterior distributions
% n -- number of samples to take

assert(numel(ts2) == numel(data2), 'inconsistent sizes between ts2 and data2');
assert(numel(ts3) == numel(data3), 'inconsistent sizes between ts3 and data3');
assert(numel(ts3) > 0, 'no suprabasal data. should run sample2');

if ~ppm_
    wbh = waitbar(0.0, 'initialising...');
end

% we use the fact that per invocation of clone_dist_* we can extract many
% separate distributions for different lambdas
mult = floor(sqrt(n)/numel(ts3));

samples_prob   = zeros(n,mult);
samples_r      = zeros(n,mult);
samples_gamma  = zeros(n,mult);
samples_lambda = zeros(n,mult);
if numel(ts2) > 0
    samples_p2     = zeros(n,mult, numel(ts2));
else
    samples_p2     = zeros(n,mult, 1);
end
samples_p3     = zeros(n,mult, numel(ts3));

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

if ~ppm_
    waitbar(0.0, wbh, 'running...');
    tstart = tic;
end

for i = 1:n
    
    [rp, r] = feval(rfun);
    lambdap = zeros(mult, 1);
    lambda  = zeros(mult, 1);
    for j = 1:mult
        [lambdap(j), lambda(j)] = feval(lambdafun);
    end
    [lambda, SI] = sort(lambda);
    lambdap = lambdap(SI);
    [gammap, gamma] = feval(gammafun);

%     fprintf(1, 'parameters: r = %f, gamma = %f\n', r, gamma);
%     fprintf(1, '            lambdas =\n');
%     fprintf(1, '                %3f\n', lambda);

    lts = zeros(numel(ts3) * mult, 1);
    tsi = zeros(numel(ts3) * mult, 1);
    for j = 1:numel(ts3)
        lts((mult * (j-1) + 1):(mult * j)) = ts3(j) .* lambda;
        tsi((mult * (j-1) + 1):(mult * j)) = j*ones(mult, 1);
    end
    [lts, SI] = sort(lts);
    tsi     = tsi(SI);
    
    dist3_ = clone_dist_bs_expv_cond(Tl, Tr, Tg, P0, r, gamma, lts);
    [m_,n_,~] = size(dist3_);
    dist3 = zeros(m_,n_,mult,numel(ts3));
    
    for j = 1:numel(ts3)
        ind = find(tsi == j);
        dist3(:,:,:,j) = dist3_(:,:,ind);
    end
    
    samples_prob(i,:)   = rp*gammap .* lambdap;
    samples_r(i,:)      = r; 
    samples_gamma(i,:)  = gamma;
    samples_lambda(i,:) = lambda; 

    if numel(ts2) > 0
        samples_p2_ = zeros(mult, numel(ts2));
        for j = 1:mult
            dist2 = clone_dist_b(r, gamma, ts2 .* lambda(j), maxN2);
            for k = 1:numel(ts2)
                % condition on b > 1
                dist2(:,k) = dist2(:,k) ./ (1 - dist2(0+1,k) - dist2(1+1,k));
                [~,col,v] = find(data2{k});
                for l = numel(v)
                    samples_p2_(j,k) = ...
                        samples_p2_(j,k) ...
                      + log(dist2(col(l),k)) * data2{k}(col(l));
                end
            end
        end
    else
        samples_p2_ = ones(mult, 1);
    end
    samples_p2(i,:,:) = samples_p2_;
    
    samples_p3_ = zeros(mult,numel(ts3));
    for j = 1:mult
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
    
    if ppm_
        ppm.increment();
    else
        waitbar(i/n, wbh, ...
            sprintf('%.1f%% complete, %.1f min to go', ...
            i/n*100, (toc(tstart)/i*n - toc(tstart))/60));
    end

end

if ~ppm_
    close(wbh);
end

samples = cell(n*mult, 6);
for i = 1:n
    for j = 1:mult
        samples{(i-1)*mult+j,1} = samples_prob(i,j);
        samples{(i-1)*mult+j,2} = samples_r(i,j);
        samples{(i-1)*mult+j,3} = samples_gamma(i,j);
        samples{(i-1)*mult+j,4} = samples_lambda(i,j);
        samples{(i-1)*mult+j,5} = squeeze(samples_p2(i,j,:));
        samples{(i-1)*mult+j,6} = squeeze(samples_p3(i,j,:));
    end
end

end
