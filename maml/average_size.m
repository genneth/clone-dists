function av = average_size(gamma, q, t)

g = sqrt(gamma * (4 * q + gamma));

av = exp(-(gamma + g) * t / 2) * ((-2 - gamma + g) + exp(g*t) * (2 + gamma + g)) / (2*g);

end