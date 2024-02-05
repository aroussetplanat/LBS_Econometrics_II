% VAR_OLS.m - VAR estimation via OLS
function results = VAR_OLS(Y, constant, p)

    [T, n] = size(Y);
    YY = Y(p+1:end,:);
    Xl = lagmatrix(Y, 1:p);

    if constant == 1
        XX = [ones(T-p, 1) Xl(p+1:end,:)];
    elseif constant == 0
        XX = [Xl(p+1:end,:)];
    else error('Error. Constant input must be 1 if a constant is used, and 0 otherwise')
    end
    
    B = XX\YY;    
    
    residuals = YY - XX*B;
    Sigma = residuals'*residuals/(T-p-n*p-1);
    results.B=B;
    results.Sigma=Sigma;
    results.residuals=residuals;

    if constant == 1
        A=B(2:end,:)';
    else
        A=B';
    end
    

    idx = [1:n:n*p, n*p+1];
    for j = 1:p
        C{j} = A(1:n, idx(j):(idx(j+1)-1));
    end
    results.C=C;

end
