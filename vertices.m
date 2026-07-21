function v = vertices(theta_min,theta_max)

% This function extracts the vertex points
% for a polytopic set of matrices

l=length(theta_min);

comb=dec2bin(0:2^l-1)' - '0';

for i=1:2^l
    for j=1:l
        if comb(j,i)==0
            v(j,i)=theta_min(j);
        else
            v(j,i)=theta_max(j);
        end
    end
end

end

