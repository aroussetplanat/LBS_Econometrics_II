% getShocks.m - Obtaining strucrural shocks
function shocks = getShocks(Y,results,A0,constant, p)

    T = size(Y,1);
    YY = Y(p+1:end,:);
    Xl = lagmatrix(Y, 1:p);
    if constant == 1
        XX = [ones(T-p, 1) Xl(p+1:end,:)];
    elseif constant == 0
        XX = [Xl(p+1:end,:)];
    else error('Error. Constant input must be 1 if a constant is used, and 0 otherwise')
    end
 
    epsilon = YY - XX*results;
   
    shocks = epsilon*inv(A0)';

end
