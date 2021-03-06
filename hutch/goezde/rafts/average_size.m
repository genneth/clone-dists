function av = average_size(r1, r2, t)

dr = r1 - r2;

av = exp(dr*t) - expm1(dr*t)*(dr-1)/dr;

av = (av - exp(-t)) / (1-exp(-t)); % condition out p_1

end
