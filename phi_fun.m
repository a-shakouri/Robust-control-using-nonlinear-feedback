function phi = phi_fun(x,z,P_total)

% Function phi(x,z) that is used in the
% construction of the feedback law

[n,l]=size(P_total);

q=l/n;

sigma=z(2);

if sigma==floor(sigma) && 1<=sigma && sigma<=q
    
    P=P_total(:,n*(sigma-1)+1:n*sigma);
    
    if x'*P*x<=z(1)
        
        phi=sigma;
        
    else
        
        phi=mod(sigma,q)+1;
        
    end
    
else
    
    phi=1;
    
end

end

