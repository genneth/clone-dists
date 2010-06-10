function condDbs3 = condDbs3(Dbs)

% drop all b=0; then condition on b+s>1

condDbs3 = Dbs;
condDbs3(0+1,:)   = 0;
condDbs3(1+1,0+1) = 0;

end
