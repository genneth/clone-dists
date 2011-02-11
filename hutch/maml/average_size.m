function av = average_size(~, gamma, ts)

% dp = gamma*d
% dd = p
% => ddp = gamma*dd = gamma*p
% => p = exp(sqrt(gamma) * t)
% => d = sqrt(gamma) * [exp(sqrt(gamma) * t) - 1]

av = exp(sqrt(gamma) * ts) + sqrt(gamma) * expm1(sqrt(gamma) * ts);

end
