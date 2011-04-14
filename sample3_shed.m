function samples = sample3_shed(rfun, gammafun, lambdafun, m, ts2, data2, ts3, data3, nsamples)

% lambdafun is unitful; ts2 and ts3 are unitful. everything else is unitless.
% rfun, gammafun, lambdafun :: (MonadRandom m) => m Double
% m :: Double -- average number of suprabasal cells per basal cell
% ts2 :: [], data2 :: {} -- times and observations of basal sizes
% ts3 :: [], data3 :: {} -- times and observations of suprabasal clone sizes
% samples :: [{r, lambda, gamma, [log p2], [log p3]}]
%     p2 and p3 are likelihoods, not posterior distributions
% nsamples -- number of samples take

assert(numel(ts2) == numel(data2), 'inconsistent sizes between ts2 and data2');
assert(numel(ts3) == numel(data3), 'inconsistent sizes between ts3 and data3');
assert(numel(ts3) > 0, 'no suprabasal data. should run sample2');

wbh = waitbar(0.0, 'initialising...', 'Name', 'sample3_shed: inference');

samples_r      = zeros(nsamples,1);
samples_gamma  = zeros(nsamples,1);
samples_lambda = zeros(nsamples,1);
if numel(ts2) > 0
    samples_p2     = zeros(nsamples,numel(ts2));
else
    samples_p2     = zeros(nsamples,1);
end
samples_p3     = zeros(nsamples,numel(ts3));

maxN2 = 0;
for i = 1:numel(ts2)
    maxN2 = max(maxN2, numel(data2{i}) - 1);
end

maxM3 = 0;
maxN3 = 0;
for i = 1:numel(ts3)
    [rows,cols,~] = find(data3{i});
    maxM3 = max(maxM3, max(rows) - 1);
    maxN3 = max(maxN3, max(cols) - 1);
end

waitbar(0.0, wbh, 'running...');
tstart = tic;

for i = 1:nsamples
    
    r = feval(rfun);
    lambda = feval(lambdafun);
    gamma = feval(gammafun);
    rho = gamma / (1 + gamma);
    mu = rho / m; % rho = mu m

    samples_r(i,:)      = r; 
    samples_gamma(i,:)  = gamma;
    samples_lambda(i,:) = lambda; 

%     fprintf(1, 'parameters: r = %f, gamma = %f\n', r, gamma);
%     fprintf(1, '            lambdas =\n');
%     fprintf(1, '                %3f\n', lambda);

    dist3 = clone_dist_bs_shed_cond(r, gamma, mu, lambda * ts3, maxM3, maxN3);
    for k = 1:numel(ts3)
        [rows,cols,v] = find(data3{k});
        for l = 1:numel(v)
            samples_p3(i,k) = samples_p3(i,k) + ...
                v(l) * log(dist3(rows(l),cols(l),k));
        end
    end
    
    % and again, with the basal data, if we have it
    if numel(ts2) > 0
        dist2 = clone_dist_b(r, gamma, lambda * ts2, maxN2);
        % condition on b>1
        [m_,~] = size(dist2);
        dist2 = dist2 ./ repmat(1 - dist2(0+1,:) - dist2(1+1,:), m_, 1);
        dist2(0+1,:) = 0; dist2(1+1,:) = 0;

        for k = 1:numel(ts2)
            [~,cols,v] = find(data2{k});
            for l = 1:numel(v)
                samples_p2(i,k) = samples_p2(i,k) + ...
                    v(l) * log(dist2(cols(l),k));
            end
        end
    else
        samples_p2(i) = 1;
    end
    
    waitbar(i/nsamples, wbh, ...
        sprintf('%.1f%% complete, %.1f min to go', ...
        i/nsamples*100, (toc(tstart)/i*nsamples - toc(tstart))/60));

end

close(wbh);

samples = cell(nsamples, 5);
for i = 1:nsamples
    samples{i,1} = samples_r(i);
    samples{i,2} = samples_gamma(i);
    samples{i,3} = samples_lambda(i);
    samples{i,4} = squeeze(samples_p2(i,:));
    samples{i,5} = squeeze(samples_p3(i,:));
end

end
