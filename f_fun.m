function f = f_fun(alpha,x,z,P_total)

% Function f for the feedback law
% which determines the state of the controller

sigma=phi_fun(x,z,P_total);

n=length(x);

P=P_total(:,n*(sigma-1)+1:n*sigma);

f=[x'*(alpha^2*P-eye(n))*x
   phi_fun(x,z,P_total)];

end

