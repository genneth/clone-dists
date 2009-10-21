function condPbs = condPbs(Pbs)

% In EYP labelling, all cells can be induced; we filter the data by
% ignoring clones of size 1 and below. The population distribution
% therefore needs to be conditioned on b+s > 1.

Z = 1 - (Pbs(0+1,0+1) + Pbs(0+1,1+1) + Pbs(1+1,0+1));
condPbs = Pbs;
condPbs(0+1,0+1) = 0;
condPbs(1+1,0+1) = 0;
condPbs(0+1,1+1) = 0;
condPbs = condPbs ./ Z;

end
