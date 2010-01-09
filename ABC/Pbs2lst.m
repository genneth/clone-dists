% converts from david's square format to the linear one
function lst = Pbs2lst(Pbs)

lst = zeros(0, 3);
[m n] = size(Pbs);

for i = 1:m
	for j = 1:n
		if Pbs(i,j) > 0
			lst(end+1,:) = [i-1 j-1 Pbs(i,j)];
		end
	end
end

end
