function condPtot = condPtot(Ptot)

Z = 1 - (Ptot(0+1) + Ptot(1+1));
condPtot = Ptot;
condPtot(0+1) = 0;
condPtot(1+1) = 0;
condPtot = condPtot ./ Z;

end
