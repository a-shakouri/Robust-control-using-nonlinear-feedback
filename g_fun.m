function g = g_fun(x,z,K_total,P_total)

% Function g for the feedback law
% which determines the control input

[~,l]=size(K_total);
n=length(x);

q=l/n;

sigma=phi_fun(x,z,P_total);

if sigma==floor(sigma) && 1<=sigma && sigma<=q
    
    g=K_total(:,n*(sigma-1)+1:n*sigma)*x;
    
else
    
    g=K_total(:,1:n)*x;
    
end

end

