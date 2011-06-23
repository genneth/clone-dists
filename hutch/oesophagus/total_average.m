function average = total_average(r, gamma, mu, lambda, ts)

raw_average = ((-1 + exp(-lambda*ts*gamma))*mu^2 - gamma*mu^2 + gamma^2*(1 - exp(-lambda*ts*mu) + mu))/(gamma*(gamma - mu)*mu);
extinction = generating_function3_shed(r, gamma, mu, lambda*ts, 0, 0, 0)';
average = raw_average ./ (1 - extinction);
% average = raw_average;

end
