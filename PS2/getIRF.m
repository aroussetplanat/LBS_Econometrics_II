% irf.m - Structural IRF Chol Identification
function IRF = getIRF(B, A0,constant, hmax, cumulateWhich)
    
     if constant == 1
        c = B(1, :)';    
        A = B(2:end,:)';
    elseif constant == 0
        A=B';
    else error('Error. Constant input must be 1 if a constant is used, and 0 otherwise')
     end

    n=size(A,1);
    p=size(A,2)/n;
    
    I = eye(n*(p-1));
    F = [A; [I zeros(n*(p-1), n)]];
    
    IRF = zeros(n, n, hmax+1);
    C = zeros(n,n,hmax+1);
    for h = 0:hmax
        Fh = F^h;
        C(:, :, h+1) = Fh(1:n, 1:n);
        IRF(:, :, h+1) = C(:, :, h+1)*A0;
    end
    IRF(cumulateWhich,:,:) = cumsum(IRF(cumulateWhich,:,:),3);
    % irf(i,j,h) shock j to var i at h
end