function [check,P,K] = qs_vertices(lambda,A_total,B_total,theta_min,theta_max)

% This function solves an LMI
% for checking quadratic stabilizability
% of a polytopic set of systems
% by only checking the vertices

l=length(theta_min);

[n,r]=size(B_total);

m=r/(l+1);

P=sdpvar(n,n);
L=sdpvar(m,n);

v=vertices(theta_min,theta_max); % Extract the vertices

Constraints=[P>=1e-6*eye(n)]; % This is to guarantee P>0

for i=1:2^l
    
    theta=v(:,i);
    
    [A,B]=sys_eval_box(A_total,B_total,theta); % Evaluate an element
                                               % within the polytpic set
                                               % corresponding to
                                               % parameter theta

    LMI=[lambda^2*P (A*P+B*L)' P
         A*P+B*L  P         zeros(n)
         P        zeros(n)  eye(n)];
    
    Constraints=[Constraints,LMI>=0];
    
end

options=sdpsettings('solver','mosek','verbose',0);

sol=optimize(Constraints,[],options);

% Check for feasibility

if length(sol.info)==27 % Checks if the problem is not feasible
                        % You may need to change 27 to another number
                        % depending on the version of your solver
    
    check=1;
    
    P=eye(n)/value(P);
    K=value(L)*P;
    
else
    
    check=0;
    
    P=zeros(n,n);
    K=zeros(m,n);
    
end



end

