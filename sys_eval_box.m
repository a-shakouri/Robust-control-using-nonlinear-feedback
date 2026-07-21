function [A,B] = sys_eval_box(A_total,B_total,theta)

% This function find an element within
% the polytoic set corresponding to
% parameter theta

[n,dn]=size(A_total);
[~,dm]=size(B_total);

l=dn/n-1;
m=dm/(l+1);

A=A_total(:,1:n);
B=B_total(:,1:m);


for i=1:l
    
    A=A+theta(i)*A_total(:,i*n+1:(i+1)*n);
    B=B+theta(i)*B_total(:,i*m+1:(i+1)*m);
    
end


end

