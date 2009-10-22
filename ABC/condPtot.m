function condPtot= condPtot(Ptot)

% In EYP labelling, all cells can be induced; we filter the data by
% ignoring clones of size 1 and below. The population distribution
% therefore needs to be conditioned on b+s > 1.

Z = 1 - (Ptot(0+1) + Ptot(1+1));
condPtot = Ptot;
condPtot(0+1) = 0;
condPtot(1+1) = 0;
condPtot = condPtot ./ Z;

end
