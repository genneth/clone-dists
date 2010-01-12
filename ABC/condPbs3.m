function condPbs3 = condPbs3(Pbs)

% drop all b=0; then condition on b+s>1

condPbs3 = Pbs;
condPbs3(0+1,:)   = 0;
condPbs3(1+1,0+1) = 0;
Z = sum(sum(condPbs3));
condPbs3 = condPbs3 ./ Z;

end
