function [c,A,C]  = reshapeVAR(beta,p,constant)
if constant == 1
    n = size(beta,2);
    c = beta(1, :)';    
    A = beta(2:end,:)'; 
    idx = [1:n:n*p, n*p+1];
    for j = 1:p
       C{j} = A(1:n, idx(j):(idx(j+1)-1));
    end
else
    n = size(beta,2);
    c = zeros(n,1)';  
    A = beta(1:end,:)';
    
    idx = [1:n:n*p, n*p+1];
    for j = 1:p
       C{j} = A(1:n, idx(j):(idx(j+1)-1));
    end
end
end