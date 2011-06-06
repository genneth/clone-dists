function samples = sample3_shed(rfun, gammafun, lambdafun, m, ts3, data3, nsamples, output_file)

% lambdafun is unitful; ts3 is unitful. everything else is unitless.
% rfun, gammafun, lambdafun :: (MonadRandom m) => m Double
% m :: Double -- average number of suprabasal cells per basal cell
% ts3 :: [], data3 :: {} -- times and observations of suprabasal clone sizes
% samples :: [{r, lambda, gamma, [log p2], [log p3]}]
%     p3 is likelihoods, not posterior distributions
% nsamples -- number of samples take

assert(numel(ts3) == numel(data3), 'inconsistent sizes between ts3 and data3');

samples = cell(nsamples, 5);

maxM3 = 0;
maxN3 = 0;
for i = 1:numel(ts3)
    [rows,cols,~] = find(data3{i});
    maxM3 = max(maxM3, max(rows) - 1);
    maxN3 = max(maxN3, max(cols) - 1);
end

start = tic;
for i = 1:nsamples
    iter_start = tic;
    
    r = feval(rfun);
    lambda = feval(lambdafun);
    gamma = feval(gammafun);
    rho = gamma / (1 + gamma);
    mu = rho / m; % rho = mu m

    samples{i,1} = r; 
    samples{i,2} = gamma;
    samples{i,3} = lambda; 

    dist3 = clone_dist_bs_shed_cond(r, gamma, mu, lambda * ts3, maxM3, maxN3);
    samples{i,5} = zeros(1,numel(ts3));
    for k = 1:numel(ts3)
        [rows,cols,v] = find(data3{k});
        for l = 1:numel(v)
            samples{i,5}(k) = samples{i,5}(k) + ...
                v(l) * log(dist3(rows(l),cols(l),k));
        end
    end
    
    % placeholder for ts2
    samples{i,4} = [];
    
    partial = samples(1:i,:);
    save(output_file, 'partial');
    
    fprintf('%d of %d, this iteration: %.1fs, to go: %.2fh\n', i, nsamples, toc(iter_start), (toc(start)/i*(nsamples-i))/3600);
end

save(output_file, 'samples');

end
